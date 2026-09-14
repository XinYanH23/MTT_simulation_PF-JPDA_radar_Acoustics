function [tracks, next_id] = track_manage(tracks, assoc, Z, R_list, ...
                                           frame_idx, next_id, cfg, sensor_types, ac_pack)
% TRACK_MANAGE  轨迹生命周期管理：更新命中/漏检计数、
%               M/N 确认、超限删除、未关联探测诞生新轨迹。
%
%   [tracks, next_id] = track_manage(tracks, assoc, Z, R_list, ...
%                                    frame_idx, next_id, cfg, sensor_types, ac_pack)
%
%   INPUTS
%     tracks    : 当前轨迹数组（含已删除轨迹）
%     assoc     : jpda_run 输出（仅活跃轨迹对应行）
%     Z         : 2×D  全部量测
%     R_list    : 2×2×D  全部量测噪声
%     frame_idx : 当前帧编号
%     next_id   : 下一条新轨迹将使用的 ID
%     cfg       : config_tracker 输出
%     sensor_types : 1×D cell（可选）
%     ac_pack   : 声学概率场（可选；用于 birth 真实性门限）
%
%   OUTPUTS
%     tracks  : 更新后轨迹数组（可能增加了新 Tentative 轨迹）
%     next_id : 更新后的全局 ID 计数器
%
%   Track Life Cycle（参考 Blackman & Popoli, §3.3）：
%
%     Tentative ──(M/N 命中)──→ Confirmed
%     Tentative ──(连续漏检 ≥ del_tentative)──→ Deleted
%     Confirmed ──(连续漏检 ≥ del_confirmed)──→ Deleted
%
%   M/N 逻辑：
%     在长度为 N 的滑动窗口 history 中，若命中次数 ≥ M，则确认。
%     此规则对抗杂波诞生的虚假 Tentative 轨迹。

D = size(Z, 2);
if nargin < 8;  sensor_types = {};  end
if nargin < 9;  ac_pack = [];       end

%% ── 阶段 1：找出活跃轨迹索引（assoc 中的行）────────────────────
if isempty(tracks) || ~isstruct(tracks)
    n_active   = 0;
    active_idx = [];
else
    active_mask = ~cellfun(@(s) strcmp(s, 'deleted'), {tracks.status});
    active_idx  = find(active_mask);
    n_active    = length(active_idx);
end

%% ── 阶段 2：更新每条活跃轨迹的命中/漏检统计 ──────────────────
for ti = 1 : n_active
    i = active_idx(ti);

    % 本帧是否有关联（β > 0 且来自至少一个门内探测）
    if isfield(assoc, 'valid_mat') && ti <= size(assoc.valid_mat, 1)
        hit = any(assoc.valid_mat(ti, :));
    else
        hit = false;
    end

    % 滑动窗口向左移一位，最右填入本帧结果
    tracks(i).history = [tracks(i).history(2:end), uint8(hit)];

    if hit
        tracks(i).consec_hits   = tracks(i).consec_hits + 1;
        tracks(i).consec_misses = 0;
        tracks(i).total_hits    = tracks(i).total_hits + 1;
        tracks(i).last_hit_frame = frame_idx;
    else
        tracks(i).consec_misses = tracks(i).consec_misses + 1;
        tracks(i).consec_hits   = 0;
        tracks(i).total_misses  = tracks(i).total_misses + 1;
    end

    tracks(i).age = tracks(i).age + 1;

    % ── 状态转移 ──────────────────────────────────────────────────
    [weak_eta, strong_eta] = ac_lifecycle_etas(cfg);
    use_del = isfield(cfg, 'ac_informed_delete') && cfg.ac_informed_delete;
    use_cfm = isfield(cfg, 'ac_aided_confirm') && cfg.ac_aided_confirm;
    La = 0;
    if (use_cfm && strcmp(tracks(i).status, 'tentative')) ...
            || (use_del && strcmp(tracks(i).status, 'confirmed') ...
                && tracks(i).consec_misses > 0)
        La = query_acoustic_likelihood_at(tracks(i).x(1:2), ac_pack, cfg);
    end

    switch tracks(i).status
        case 'tentative'
            hits_in_window = sum(tracks(i).history);
            M_need = cfg.mn_M;
            if use_cfm && La >= strong_eta
                M_need = cfg.mn_M_if_strong;
            end
            if hits_in_window >= M_need
                tracks(i).status = 'confirmed';
                if cfg.debug
                    fprintf('[帧%4d] 轨迹 #%d 确认 (命中 %d/%d, L_a=%.2e)\n', ...
                        frame_idx, tracks(i).id, hits_in_window, cfg.mn_N, La);
                end
            end
            if tracks(i).consec_misses >= cfg.del_tentative
                tracks(i).status = 'deleted';
                if cfg.debug
                    fprintf('[帧%4d] 轨迹 #%d Tentative 删除\n', ...
                        frame_idx, tracks(i).id);
                end
            end

        case 'confirmed'
            del_now = tracks(i).consec_misses >= cfg.del_confirmed;
            if ~del_now && use_del
                early = cfg.del_confirmed_if_weak;
                if tracks(i).consec_misses >= early && La < weak_eta
                    del_now = true;
                end
            end
            if del_now
                tracks(i).status = 'deleted';
                if cfg.debug
                    fprintf('[帧%4d] 轨迹 #%d Confirmed 删除 (连续漏检%d帧, L_a=%.2e)\n', ...
                        frame_idx, tracks(i).id, tracks(i).consec_misses, La);
                end
            end
    end

    % 记录历史
    tracks(i).state_hist(:, end+1)   = tracks(i).x;
    tracks(i).P_hist(:,:,end+1)      = tracks(i).P;
    tracks(i).mu_hist(:, end+1)      = tracks(i).mu;
    tracks(i).time_hist(end+1)       = frame_idx;
end

%% ── 阶段 3：找出完全未关联的探测，诞生新轨迹 ─────────────────
% "完全未关联"：该探测不在任何活跃轨迹的门内
if n_active > 0 && isfield(assoc, 'valid_mat') && size(assoc.valid_mat, 2) == D
    unassigned_dets = find(~any(assoc.valid_mat, 1));
else
    unassigned_dets = 1 : D;    % 没有活跃轨迹，所有探测均未关联
end

allow_birth = (~isfield(cfg, 'allow_birth')) || cfg.allow_birth;
radar_only  = isfield(cfg, 'birth_radar_only') && cfg.birth_radar_only;
suppress_r  = 0;
if isfield(cfg, 'birth_suppress_radius') && ~isempty(cfg.birth_suppress_radius)
    suppress_r = cfg.birth_suppress_radius;
end
eta_a = 0;
if isfield(cfg, 'birth_acoustic_eta') && ~isempty(cfg.birth_acoustic_eta)
    eta_a = cfg.birth_acoustic_eta;
end

% 仅对已确认航迹做邻近抑制（tentative 杂波不应封锁真实新目标）
conf_pos = zeros(2, 0);
if n_active > 0
    for ti = 1 : n_active
        i = active_idx(ti);
        if strcmp(tracks(i).status, 'confirmed')
            conf_pos = [conf_pos, tracks(i).x(1:2)]; %#ok<AGROW>
        end
    end
end

% 本帧已接受的 birth 位置（同帧去重）
born_pos = zeros(2, 0);

for d = unassigned_dets
    if ~allow_birth
        continue
    end
    % 仅雷达诞生：跳过非雷达（声学/未知）未关联探测，避免冗余轨迹爆发
    if radar_only && ~isempty(sensor_types) && d <= numel(sensor_types) ...
            && ~strcmpi(sensor_types{d}, 'Radar')
        continue
    end
    z = Z(:, d);
    R = R_list(:, :, d);

    % 邻近抑制：落在已确认航迹附近的未关联点不新建航迹
    if suppress_r > 0 && ~isempty(conf_pos)
        dmin = min(sqrt(sum((conf_pos - z).^2, 1)));
        if dmin <= suppress_r
            if cfg.debug
                fprintf('[帧%4d] 抑制诞生 at (%.1f, %.1f)，距最近确认航迹 %.1f m\n', ...
                    frame_idx, z(1), z(2), dmin);
            end
            continue
        end
    end
    % 同帧 birth 去重
    if suppress_r > 0 && ~isempty(born_pos)
        dmin = min(sqrt(sum((born_pos - z).^2, 1)));
        if dmin <= suppress_r
            continue
        end
    end

    % 声学真实性门限（论文 3.2.4）：L_a(z) > eta_a
    if eta_a > 0 && ~isempty(ac_pack) && isstruct(ac_pack) ...
            && isfield(ac_pack, 'detected') && ac_pack.detected
        La = query_acoustic_likelihood_at(z, ac_pack, cfg);
        if La < eta_a
            if cfg.debug
                fprintf('[帧%4d] 声学门限拒绝诞生 at (%.1f, %.1f), L_a=%.2e < eta=%.2e\n', ...
                    frame_idx, z(1), z(2), La, eta_a);
            end
            continue
        end
    end

    new_track = track_create(z, R, next_id, frame_idx, cfg);
    next_id   = next_id + 1;
    born_pos  = [born_pos, z]; %#ok<AGROW>

    % 追加到轨迹数组（对齐字段：主循环帧末会 rmfield c_bar_，
    % 而 track_create 仍带该字段，直接追加会触发“结构体字段不一致”）
    if isempty(tracks)
        tracks = new_track;
    else
        new_track = align_track_struct(new_track, tracks(1));
        tracks(end + 1) = new_track;  %#ok<AGROW>
    end

    if cfg.debug
        fprintf('[帧%4d] 诞生新轨迹 #%d at (%.1f, %.1f)\n', ...
            frame_idx, new_track.id, z(1), z(2));
    end
end
end

% -------------------------------------------------------------------------
function t = align_track_struct(t, ref)
%ALIGN_TRACK_STRUCT  使 t 与 ref 字段集合一致，便于 struct 数组追加。
fn_t   = fieldnames(t);
fn_ref = fieldnames(ref);
for i = 1:numel(fn_ref)
    f = fn_ref{i};
    if ~isfield(t, f)
        t.(f) = ref.(f);   % 占位；随后会被正常流程覆盖
        if isnumeric(t.(f)); t.(f) = []; end
        if iscell(t.(f));    t.(f) = {}; end
    end
end
extra = setdiff(fn_t, fn_ref);
if ~isempty(extra)
    t = rmfield(t, extra);
end
end

function [weak, strong] = ac_lifecycle_etas(cfg)
weak = 1e-3;
if isfield(cfg, 'birth_acoustic_eta') && ~isempty(cfg.birth_acoustic_eta)
    weak = cfg.birth_acoustic_eta;
end
if isfield(cfg, 'ac_field_weak_eta') && ~isempty(cfg.ac_field_weak_eta)
    weak = cfg.ac_field_weak_eta;
end
if isfield(cfg, 'ac_field_strong_eta') && ~isempty(cfg.ac_field_strong_eta)
    strong = cfg.ac_field_strong_eta;
else
    strong = 1e-2;   % Softmax 真值处约 2e-2
end
end

function La = query_acoustic_likelihood_at(z, ac_pack, cfg)
%QUERY_ACOUSTIC_LIKELIHOOD_AT  在位置 z=[px;py](m) 查询声学似然场。
if nargin < 3 || isempty(cfg); cfg = struct(); end
if isempty(ac_pack) || ~isstruct(ac_pack) ...
        || ~isfield(ac_pack, 'detected') || ~ac_pack.detected
    La = 1e-6;
    return
end
kappa = 1000;
if isfield(cfg, 'pf_m2ac_scale') && ~isempty(cfg.pf_m2ac_scale)
    kappa = cfg.pf_m2ac_scale;
end
bg = 1e-6;
if isfield(cfg, 'pf_ac_bg_likelihood') && ~isempty(cfg.pf_ac_bg_likelihood)
    bg = cfg.pf_ac_bg_likelihood;
end

% 生命周期/诞生门限一律查 Softmax（与 η_a 同一标度）。L_rel 只用于 PF 权重。
L_field = [];
method_ac = 'softmax';
if isfield(ac_pack, 'method') && ~isempty(ac_pack.method)
    method_ac = lower(ac_pack.method);
end
switch method_ac
    case 'intensity'
        if isfield(ac_pack, 'L_intensity'); L_field = ac_pack.L_intensity; end
    case 'nms_gmm'
        if isfield(ac_pack, 'L_nms_gmm'); L_field = ac_pack.L_nms_gmm; end
    otherwise
        if isfield(ac_pack, 'L_softmax'); L_field = ac_pack.L_softmax; end
end
if isempty(L_field) || ~isfield(ac_pack, 'xc') || ~isfield(ac_pack, 'zc')
    La = bg;
    return
end

px_mm = z(1) * kappa;
py_mm = z(2) * kappa;
try
    F = griddedInterpolant({ac_pack.xc(:), ac_pack.zc(:)}, double(L_field), 'linear', 'none');
    La = F(px_mm, py_mm);
    if isnan(La); La = bg; end
catch
    La = bg;
end
La = max(La, bg);
end
