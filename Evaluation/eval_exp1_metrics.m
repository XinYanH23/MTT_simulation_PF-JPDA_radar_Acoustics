function m = eval_exp1_metrics(tracks, truth, time, cfg, frame_log, quiet)
%EVAL_EXP1_METRICS  实验一专用跟踪级指标。
%
%   重点：虚警抑制、基数估计、航迹连续性（非定位精度）。
%   quiet=true 时关闭逐次 fprintf（蒙特卡洛批量用）。

if nargin < 6 || isempty(quiet)
    quiet = false;
end

assoc_radius = 50;   % m
c_ospa = 100;
p_ospa = 2;
N = numel(time);
Nt = numel(truth);

% 曾确认航迹
conf_idx = [];
for i = 1:numel(tracks)
    if strcmp(tracks(i).status, 'confirmed') || tracks(i).total_hits >= cfg.mn_M
        conf_idx(end+1) = i; %#ok<AGROW>
    end
end
nc = numel(conf_idx);

%% ── 逐帧：真值/估计快照 + 航迹级关联（供碎裂 / ID / OSPA 分解）──
ospa_k     = nan(1, N);
ospa_loc_k = nan(1, N);
ospa_card_k = nan(1, N);
ospa_assign_k = nan(1, N);
card_err_k = nan(1, N);
n_true_k   = zeros(1, N);
n_est_k    = zeros(1, N);
frame_assign = nan(Nt, N);      % 真目标 j 在帧 k 匹配的航迹 ID（50 m）
track_hits_truth = false(nc, Nt);
gt_xy  = cell(1, N);
gt_ids = cell(1, N);
est_xy = cell(1, N);
est_ids_k = cell(1, N);

id2ci = containers.Map('KeyType', 'double', 'ValueType', 'double');
for ci = 1:nc
    id2ci(tracks(conf_idx(ci)).id) = ci;
end

for k = 1:N
    gt = zeros(2, 0);
    gt_id = [];
    for j = 1:Nt
        raw = truth(j).pos3D;
        if k <= size(raw,1) && ~any(isnan(raw(k,1:2)))
            gt(:, end+1) = raw(k,1:2)'; %#ok<AGROW>
            gt_id(end+1) = j; %#ok<AGROW>
        end
    end
    n_true_k(k) = size(gt, 2);

    if isfield(frame_log, 'n_confirmed')
        n_est_k(k) = frame_log.n_confirmed(k);
    end
    card_err_k(k) = abs(n_est_k(k) - n_true_k(k));

    if isfield(frame_log, 'conf_pos') && ~isempty(frame_log.conf_pos{k})
        est = frame_log.conf_pos{k};
        est_ids = frame_log.conf_ids{k};
    else
        est = zeros(2, 0);
        est_ids = [];
    end
    gt_xy{k} = gt;  gt_ids{k} = gt_id;
    est_xy{k} = est; est_ids_k{k} = est_ids;

    if n_true_k(k) > 0 && ~isempty(est)
        D = pdist2(gt', est');
        Dcut = D;
        Dcut(D > assoc_radius) = 1e6;
        [r, c] = match_rows(Dcut);
        for t = 1:numel(r)
            if Dcut(r(t), c(t)) < 1e5
                j = gt_id(r(t));
                tid = est_ids(c(t));
                frame_assign(j, k) = tid;
                if isKey(id2ci, tid)
                    track_hits_truth(id2ci(tid), j) = true;
                end
            end
        end
    end
end

% 每个真目标的主航迹 = 关联帧上出现次数最多的确认 ID（与碎裂/IDS 同源）
primary_id = nan(Nt, 1);
for j = 1:Nt
    seq = frame_assign(j, :);
    seq = seq(~isnan(seq));
    if ~isempty(seq)
        primary_id(j) = mode(seq);
    end
end

for k = 1:N
    [ospa_k(k), ospa_loc_k(k), ospa_card_k(k), ospa_assign_k(k)] = ...
        ospa_decompose(gt_xy{k}, est_xy{k}, gt_ids{k}, est_ids_k{k}, ...
                       primary_id, c_ospa, p_ospa);
end

active = n_true_k > 0;
valid_ospa = active & ~isnan(ospa_k);
m.ospa_k        = ospa_k;
m.ospa_mean     = mean(ospa_k(valid_ospa));
m.ospa_rmse     = sqrt(mean(ospa_k(valid_ospa).^2));
m.ospa_loc_k    = ospa_loc_k;
m.ospa_card_k   = ospa_card_k;
m.ospa_assign_k = ospa_assign_k;
m.ospa_loc_mean    = mean(ospa_loc_k(valid_ospa));
m.ospa_card_mean   = mean(ospa_card_k(valid_ospa));
m.ospa_assign_mean = mean(ospa_assign_k(valid_ospa));
m.card_err_k    = card_err_k;
m.card_err_mean = mean(card_err_k(active));
m.n_true_k      = n_true_k;
m.n_est_k       = n_est_k;
m.peak_confirmed = max(n_est_k);

%% ── 虚警确认比例 ─────────────────────────────────────────────────
clutter_mask = ~any(track_hits_truth, 2);
n_false = sum(clutter_mask);
m.n_confirmed_total   = nc;
m.n_false_confirmed   = n_false;
m.false_confirm_ratio = n_false / max(nc, 1);

%% ── 航迹寿命 ─────────────────────────────────────────────────────
lives = zeros(1, nc);
for ci = 1:nc
    i = conf_idx(ci);
    lives(ci) = tracks(i).time_hist(end) - tracks(i).birth_frame + 1;
end
m.track_lives       = lives;
m.mean_track_life   = mean(lives);
m.median_track_life = median(lives);

%% ── 碎裂 ─────────────────────────────────────────────────────────
frag = sum(track_hits_truth, 1);
m.fragmentation_per_truth = frag;
m.mean_frag_per_truth     = mean(frag);
m.fragmentation_count     = sum(max(frag - 1, 0));

%% ── ID Switch ────────────────────────────────────────────────────
id_sw = 0;
for j = 1:Nt
    prev = NaN;
    for k = 1:N
        cur = frame_assign(j, k);
        if isnan(cur), continue; end
        if ~isnan(prev) && cur ~= prev
            id_sw = id_sw + 1;
        end
        prev = cur;
    end
end
m.id_switches = id_sw;

%% ── 附带位置 RMSE ────────────────────────────────────────────────
base = eval_metrics(tracks, truth, time, cfg, quiet);
m.pos_rmse = base.pos_rmse;

if ~quiet
    fprintf('\n===== 实验一指标 =====\n');
    fprintf('  OSPA mean / RMSE     : %.2f / %.2f m\n', m.ospa_mean, m.ospa_rmse);
    fprintf('  OSPA loc/card/assign : %.2f / %.2f / %.2f m\n', ...
        m.ospa_loc_mean, m.ospa_card_mean, m.ospa_assign_mean);
    fprintf('  基数误差均值         : %.2f\n', m.card_err_mean);
    fprintf('  峰值确认航迹数       : %d\n', m.peak_confirmed);
    fprintf('  虚警确认比例         : %.1f%% (%d/%d)\n', ...
        100*m.false_confirm_ratio, m.n_false_confirmed, m.n_confirmed_total);
    fprintf('  平均/中位航迹寿命    : %.1f / %.1f 帧\n', m.mean_track_life, m.median_track_life);
    fprintf('  碎裂次数 / 每目标均值: %d / %.2f\n', m.fragmentation_count, m.mean_frag_per_truth);
    fprintf('  ID Switch            : %d\n', m.id_switches);
    fprintf('  位置 RMSE（附带）    : %.2f m\n', m.pos_rmse);
    fprintf('======================\n\n');
end
end

function [d, d_loc, d_card, d_asg] = ospa_decompose(X, Y, x_id, y_id, primary_id, c, p)
%OSPA_DECOMPOSE  将 OSPA 拆成定位 / 基数 / 分配三项（p 范数可加）。
%
%   匈牙利配对（截断距离 c）与标准 OSPA 相同，因此
%     OSPA^p = OSPA_loc^p + OSPA_card^p + OSPA_assign^p。
%   身份一致（估计 ID = 该真目标主航迹）的配对计入 loc；
%   身份不一致的配对计入 assign（对应 ID Switch 与航迹碎裂）；
%   |nX-nY| 个截断虚拟配对计入 card。

n = size(X, 2); m = size(Y, 2);
nn = max(n, m);
if n == 0 && m == 0
    d = 0; d_loc = 0; d_card = 0; d_asg = 0;
    return
end
if n == 0 || m == 0
    d = c; d_loc = 0; d_card = c; d_asg = 0;
    return
end

D = min(pdist2(X', Y'), c);
pairs = match_pairs(D.^p);
sum_loc = 0;
sum_asg = 0;
for t = 1:size(pairs, 1)
    di  = D(pairs(t,1), pairs(t,2));
    gj  = x_id(pairs(t,1));
    tid = y_id(pairs(t,2));
    pid = primary_id(gj);
    if ~isnan(pid) && tid ~= pid
        sum_asg = sum_asg + di^p;
    else
        sum_loc = sum_loc + di^p;
    end
end
sum_card = c^p * abs(n - m);
d      = ((sum_loc + sum_asg + sum_card) / nn)^(1/p);
d_loc  = (sum_loc  / nn)^(1/p);
d_card = (sum_card / nn)^(1/p);
d_asg  = (sum_asg  / nn)^(1/p);
end

function pairs = match_pairs(Cost)
% 返回 P×2 配对 [row, col]
try
    pairs = matchpairs(Cost, 1e7);
catch
    C = Cost;
    pairs = zeros(0, 2);
    while true
        [v, idx] = min(C(:));
        if ~isfinite(v) || v >= 1e6, break; end
        [ri, ci] = ind2sub(size(C), idx);
        pairs(end+1, :) = [ri, ci]; %#ok<AGROW>
        C(ri, :) = Inf; C(:, ci) = Inf;
        if size(pairs,1) >= min(size(Cost)), break; end
    end
end
end

function [r, c] = match_rows(Cost)
pairs = match_pairs(Cost);
if isempty(pairs)
    r = []; c = [];
else
    r = pairs(:,1); c = pairs(:,2);
end
end
