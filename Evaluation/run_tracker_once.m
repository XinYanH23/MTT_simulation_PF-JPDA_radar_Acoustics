function out = run_tracker_once(allDetections, truth, time, cfg, ac_cfg, opts)
%RUN_TRACKER_ONCE  单次 IMM-JPDA-(PF|KF) 跟踪（供对比实验复用）。
%
%   opts.enable_acoustic : 兼容旧接口；true→field，false→none
%   opts.acoustic_mode   : 'none' | 'field' | 'point_map' | 'point_multi'
%       none        — Radar-only
%       field       — 声学概率场：KF/PF 共用同一 ac_pack 做更新与 η_a 起始
%       point_map   — 场退化为单点（MAP）进入 JPDA，更新不用场
%       point_multi — 场多峰（局部极大）进入 JPDA，更新不用场
%   opts.verbose, opts.method, opts.point_K（多峰上限，默认 3）
%   cfg.use_pf           : true=IMM-PF；false=IMM-KF（须由调用方设置）

if nargin < 6 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'verbose'); opts.verbose = true; end
if ~isfield(opts, 'method');  opts.method = 'softmax'; end
if ~isfield(opts, 'point_K'); opts.point_K = 3; end

if ~isfield(opts, 'acoustic_mode') || isempty(opts.acoustic_mode)
    if isfield(opts, 'enable_acoustic') && ~opts.enable_acoustic
        opts.acoustic_mode = 'none';
    else
        opts.acoustic_mode = 'field';
    end
end
mode = lower(string(opts.acoustic_mode));

N = length(time);
cfg.dt = time(2) - time(1);
if ~isfield(cfg, 'use_pf') || isempty(cfg.use_pf)
    cfg.use_pf = true;
end
use_pf = logical(cfg.use_pf);

% 点量测模式：关闭声学起始门限（改由点迹进 JPDA）；场模式保留 η_a
if mode == "point_map" || mode == "point_multi" || mode == "none"
    cfg.birth_acoustic_eta = 0;
end

tracks   = struct([]);
next_id  = 1;
prev_post_ac = [];
ac_pack_hist = cell(N, 1);
assoc_hist   = cell(N, 1);
frame_log = struct();
frame_log.n_confirmed = zeros(1, N);
frame_log.n_tentative = zeros(1, N);
frame_log.n_deleted   = zeros(1, N);
frame_log.conf_pos    = cell(1, N);
frame_log.conf_ids    = cell(1, N);

if opts.verbose
    filt = ternary_str(use_pf, 'PF', 'KF');
    fprintf('===== 跟踪开始 [%s | %s] 帧数=%d =====\n', mode, filt, N);
end

for k = 1:N
    if opts.verbose && mod(k, 50) == 0
        if isempty(tracks)
            n_act = 0;
        else
            n_act = sum(~strcmp({tracks.status}, 'deleted'));
        end
        fprintf('  帧 %4d/%d | 活跃: %d\n', k, N, n_act);
    end

    %% Step 0：声学场（none 跳过；field/point 均生成）
    ac_pack_k = [];
    need_field = (mode ~= "none") && ~isempty(ac_cfg);
    if need_field
        src_pos_k = [];
        for tn = 1:length(truth)
            if k <= size(truth(tn).pos3D, 1) && ~any(isnan(truth(tn).pos3D(k,:)))
                src_pos_k = [src_pos_k; truth(tn).pos3D(k, 1:3)]; %#ok<AGROW>
            end
        end
        if ~isempty(src_pos_k)
            [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
                src_pos_k, ac_cfg, prev_post_ac, opts.method);
            ac_pack_hist{k} = ac_pack_k;
        end
    end

    %% Step 1：IMM 混合 + 预测（PF / KF）
    for i = 1:length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end
        if use_pf
            [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
                                                       tracks(i).mu, tracks(i).Pi);
            tracks(i).c_bar_ = c_bar_i;
            tracks(i).models = imm_pf_predict(tracks(i).models);
        else
            [tracks(i).models, c_bar_i] = imm_mix(tracks(i).models, ...
                                                   tracks(i).mu, tracks(i).Pi);
            tracks(i).c_bar_ = c_bar_i;
            tracks(i).models = imm_predict(tracks(i).models);
        end
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% Step 2&3：量测 + JPDA
    [Z_k, R_list_k, sensor_types_k] = detections_to_ZR(allDetections{k}, cfg);
    if (~isfield(cfg,'jpda_radar_only')) || cfg.jpda_radar_only
        if ~isempty(sensor_types_k)
            is_radar       = cellfun(@(s) strcmpi(s,'Radar'), sensor_types_k);
            Z_k            = Z_k(:, is_radar);
            R_list_k       = R_list_k(:, :, is_radar);
            sensor_types_k = sensor_types_k(is_radar);
        end
    end

    % 点量测模式：把声学场峰值并入 JPDA 量测集
    if (mode == "point_map" || mode == "point_multi") && ~isempty(ac_pack_k) ...
            && isstruct(ac_pack_k) && isfield(ac_pack_k,'detected') && ac_pack_k.detected
        [Z_ac, R_ac, types_ac] = acoustic_field_to_points(ac_pack_k, cfg, mode, opts.point_K);
        if ~isempty(Z_ac)
            Z_k = [Z_k, Z_ac];
            R_list_k = cat(3, R_list_k, R_ac);
            sensor_types_k = [sensor_types_k, types_ac];
        end
    end

    if isempty(tracks)
        active_mask = false(1, 0);
    else
        active_mask = ~strcmp({tracks.status}, 'deleted');
    end
    active_idx = find(active_mask);
    n_active   = length(active_idx);

    if n_active > 0 && size(Z_k, 2) > 0
        assoc = jpda_run(tracks(active_idx), Z_k, R_list_k, cfg);
    else
        assoc.valid_mat  = false(n_active, size(Z_k,2));
        assoc.innov_data = cell(n_active, size(Z_k,2));
        assoc.beta       = zeros(n_active, size(Z_k,2));
        assoc.beta0      = ones(n_active, 1);
        assoc.clusters   = {};
    end
    assoc_hist{k} = assoc;

    %% Step 4：PDA 更新（field 模式下 KF/PF 均传入同一 ac_pack）
    ac_for_update = [];
    if mode == "field"
        ac_for_update = ac_pack_k;
    end
    for ti = 1:n_active
        i = active_idx(ti);
        tracks(i) = track_imm_pda_update(tracks(i), Z_k, R_list_k, ...
                                          assoc, ti, cfg, ac_for_update);
    end

    %% Step 5：轨迹管理（field 模式：KF/PF 均可用 η_a 起始门限）
    ac_for_birth = [];
    if mode == "field"
        ac_for_birth = ac_pack_k;
    end
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, ...
                                      k, next_id, cfg, sensor_types_k, ac_for_birth);
    tracks = track_merge(tracks, cfg);

    if ~isempty(tracks)
        st = {tracks.status};
        frame_log.n_confirmed(k) = sum(strcmp(st,'confirmed'));
        frame_log.n_tentative(k) = sum(strcmp(st,'tentative'));
        frame_log.n_deleted(k)   = sum(strcmp(st,'deleted'));
        conf_i = find(strcmp(st, 'confirmed'));
        if ~isempty(conf_i)
            pos = zeros(2, numel(conf_i));
            ids = zeros(1, numel(conf_i));
            for ii = 1:numel(conf_i)
                pos(:, ii) = tracks(conf_i(ii)).x(1:2);
                ids(ii)    = tracks(conf_i(ii)).id;
            end
            frame_log.conf_pos{k} = pos;
            frame_log.conf_ids{k} = ids;
        else
            frame_log.conf_pos{k} = zeros(2, 0);
            frame_log.conf_ids{k} = [];
        end
    else
        frame_log.conf_pos{k} = zeros(2, 0);
        frame_log.conf_ids{k} = [];
    end
    if ~isempty(tracks) && isfield(tracks,'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end
end

if opts.verbose
    fprintf('===== 跟踪完成 =====\n');
end

metrics = eval_metrics(tracks, truth, time, cfg);

out = struct();
out.tracks       = tracks;
out.metrics      = metrics;
out.frame_log    = frame_log;
out.ac_pack_hist = ac_pack_hist;
out.assoc_hist   = assoc_hist;
out.cfg          = cfg;
out.acoustic_mode = char(mode);
out.enable_acoustic = (mode == "field");
out.use_pf = use_pf;
end

function s = ternary_str(cond, a, b)
if cond, s = a; else, s = b; end
end

% -------------------------------------------------------------------------
function [Z_ac, R_ac, types_ac] = acoustic_field_to_points(ac_pack, cfg, mode, K)
Z_ac = zeros(2, 0);
R_ac = zeros(2, 2, 0);
types_ac = {};

if ~isfield(ac_pack, 'L_softmax') || isempty(ac_pack.L_softmax)
    return
end
L = ac_pack.L_softmax;
if isfield(ac_pack, 'XX') && isfield(ac_pack, 'ZZ') && isequal(size(ac_pack.XX), size(L))
    XX = ac_pack.XX; ZZ = ac_pack.ZZ;
else
    xc = ac_pack.xc(:); zc = ac_pack.zc(:);
    [XX, ZZ] = meshgrid(xc, zc);
    if ~isequal(size(L), size(XX))
        [XX, ZZ] = ndgrid(xc, zc);
    end
    if ~isequal(size(L), size(XX))
        return
    end
end

Rdef = cfg.R_acoustic;
if mode == "point_map"
    [~, est_map] = field_point_estimate(L, XX, ZZ);  % mm
    Z_ac = [est_map(1); est_map(2)] / 1000;
    R_ac = Rdef;
    types_ac = {'Acoustic'};
    return
end

peaks = local_maxima_2d(L, K);
for i = 1:size(peaks, 1)
    r = peaks(i,1); c = peaks(i,2);
    z = [XX(r, c); ZZ(r, c)] / 1000;
    if any(~isfinite(z)), continue; end
    Z_ac(:, end+1) = z; %#ok<AGROW>
    R_ac(:, :, end+1) = Rdef; %#ok<AGROW>
    types_ac{end+1} = 'Acoustic'; %#ok<AGROW>
end

if isempty(Z_ac)
    [~, est_map] = field_point_estimate(L, XX, ZZ);
    Z_ac = [est_map(1); est_map(2)] / 1000;
    R_ac = Rdef;
    types_ac = {'Acoustic'};
end
end

function peaks = local_maxima_2d(L, K)
[nr, nc] = size(L);
mask = false(nr, nc);
for i = 2:nr-1
    for j = 2:nc-1
        v = L(i,j);
        if v >= L(i-1,j) && v >= L(i+1,j) && v >= L(i,j-1) && v >= L(i,j+1) ...
                && v >= L(i-1,j-1) && v >= L(i-1,j+1) && v >= L(i+1,j-1) && v >= L(i+1,j+1)
            mask(i,j) = true;
        end
    end
end
[ii, jj] = find(mask);
if isempty(ii)
    [~, idx] = max(L(:));
    [r, c] = ind2sub(size(L), idx);
    peaks = [r, c];
    return
end
vals = L(sub2ind(size(L), ii, jj));
[~, ord] = sort(vals, 'descend');
n = min(K, numel(ord));
peaks = [ii(ord(1:n)), jj(ord(1:n))];
end
