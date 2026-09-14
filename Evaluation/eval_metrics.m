function metrics = eval_metrics(tracks, truth, time, cfg, quiet)
% EVAL_METRICS  计算多目标跟踪性能指标。
%
%   metrics = eval_metrics(tracks, truth, time, cfg)
%   metrics = eval_metrics(tracks, truth, time, cfg, quiet)
%
%   INPUTS
%     tracks  : 跟踪器输出的轨迹数组（含 state_hist, time_hist, status）
%     truth   : Step1/Step2 真值结构体数组（含 pos3D）
%     time    : 1×N 时间向量
%     cfg     : config_tracker 输出
%
%   OUTPUT
%     metrics 结构体，含以下字段：
%
%     .pos_rmse         平均位置 RMSE (m)
%     .vel_rmse         平均速度 RMSE (m/s)
%     .confirmed_tracks 最终确认轨迹数
%     .false_tracks     确认后又删除的虚假轨迹数（近似）
%     .confirm_delay    平均确认延迟（帧）
%     .track_frag       轨迹碎片化率（已确认轨迹中途删除并重建的比例）
%     .association_rate 活跃帧的平均关联率
%
%   评估方法：
%     使用贪心最小距离赋值将估计轨迹与真值轨迹对应，
%     在各帧上计算 RMSE。仅统计目标真实存在（pos3D 非 NaN）的帧，
%     不对异步出生/死亡区间的真值做插值或外推。

if nargin < 5 || isempty(quiet)
    quiet = false;
end

N_truth = length(truth);
N_frame = length(time);

%% ── 筛选已确认（或曾经确认）轨迹 ──────────────────────────────
conf_mask = false(1, length(tracks));
for i = 1 : length(tracks)
    if strcmp(tracks(i).status, 'confirmed') || tracks(i).total_hits >= cfg.mn_M
        conf_mask(i) = true;
    end
end
conf_idx = find(conf_mask);

metrics.confirmed_tracks = sum(conf_mask);
metrics.false_tracks     = 0;   % 见下
metrics.confirm_delay    = NaN;
metrics.pos_rmse         = NaN;
metrics.vel_rmse         = NaN;
metrics.track_frag       = NaN;
metrics.association_rate = NaN;
metrics.rmse_pos_k       = nan(1, N_frame);
metrics.rmse_vel_k       = nan(1, N_frame);
metrics.ospa             = nan(1, N_frame);

if isempty(conf_idx) || N_truth == 0
    fprintf('[eval_metrics] 无已确认轨迹或无真值，跳过 RMSE 计算。\n');
    return
end

%% ── 平均确认延迟 ────────────────────────────────────────────────
delays = [];
for i = conf_idx
    if ~isnan(tracks(i).birth_frame)
        delays(end+1) = tracks(i).last_hit_frame - tracks(i).birth_frame + 1; %#ok<AGROW>
    end
end
if ~isempty(delays)
    metrics.confirm_delay = mean(delays);
end

%% ── 位置 / 速度 RMSE（逐帧 + 全局均值）────────────────────────────
% 贪心赋值：对每帧，用最近邻将已确认轨迹与真值对应
pos_err_sq  = [];
c_ospa      = 100;
p_ospa      = 2;

for k = 1 : N_frame
    % 获取本帧各已确认轨迹的估计位置
    est_pos = zeros(2, length(conf_idx));
    est_vel = zeros(2, length(conf_idx));
    valid_est = false(1, length(conf_idx));

    for ci = 1 : length(conf_idx)
        i = conf_idx(ci);
        t_hist = tracks(i).time_hist;
        [~, loc] = min(abs(t_hist - k));
        if ~isempty(loc) && abs(t_hist(loc) - k) <= 1
            s = tracks(i).state_hist(:, loc);
            est_pos(:, ci)  = s(1:2);
            est_vel(:, ci)  = s(3:4);
            valid_est(ci)   = true;
        end
    end

    % 获取本帧活跃真值（不对 NaN 帧外推）
    gt_pos = zeros(2, 0);
    for j = 1 : N_truth
        raw = truth(j).pos3D;
        if k <= size(raw, 1) && ~any(isnan(raw(k, 1:2)))
            gt_pos(:, end+1) = raw(k, 1:2)'; %#ok<AGROW>
        end
    end
    if isempty(gt_pos)
        continue
    end

    % 贪心最近邻赋值（仅活跃真值）
    est_used = false(1, length(conf_idx));
    n_active = size(gt_pos, 2);
    frame_err_sq = [];

    for j = 1 : n_active
        best_ci  = -1;
        best_d   = inf;
        for ci = 1 : length(conf_idx)
            if ~valid_est(ci) || est_used(ci), continue; end
            d = norm(est_pos(:,ci) - gt_pos(:,j));
            if d < best_d
                best_d  = d;
                best_ci = ci;
            end
        end
        if best_ci > 0 && best_d < 100   % 100 m 距离上限
            pos_err_sq(end+1) = best_d^2;      %#ok<AGROW>
            frame_err_sq(end+1) = best_d^2;    %#ok<AGROW>
            est_used(best_ci) = true;
        end
    end

    if ~isempty(frame_err_sq)
        metrics.rmse_pos_k(k) = sqrt(mean(frame_err_sq));
    end

    % 逐帧 OSPA（p=2, c=100m）
    est_xy = est_pos(:, valid_est);
    if isempty(gt_pos)
        metrics.ospa(k) = c_ospa;
    elseif isempty(est_xy)
        metrics.ospa(k) = c_ospa;
    else
        m = size(est_xy, 2);
        n_t = size(gt_pos, 2);
        mn = max(m, n_t);
        D = zeros(m, n_t);
        for ii = 1:m
            for jj = 1:n_t
                D(ii,jj) = min(norm(est_xy(:,ii) - gt_pos(:,jj)), c_ospa);
            end
        end
        asgn = 0;
        D2 = D;
        for s = 1:min(m, n_t)
            [~, bi] = min(D2(:));
            [r, c_] = ind2sub(size(D2), bi);
            asgn = asgn + D2(r, c_)^p_ospa;
            D2(r,:) = Inf; D2(:, c_) = Inf;
        end
        metrics.ospa(k) = ((asgn + c_ospa^p_ospa * abs(m - n_t)) / mn)^(1/p_ospa);
    end
end

if ~isempty(pos_err_sq)
    metrics.pos_rmse = sqrt(mean(pos_err_sq));
end
metrics.rmse_pos = metrics.pos_rmse;

%% ── 虚假轨迹率 ─────────────────────────────────────────────────
% 确认后被删除、且存在时间短（< N 帧）的视为虚假轨迹
for i = conf_idx
    if strcmp(tracks(i).status, 'deleted') && tracks(i).age < cfg.mn_N * 3
        metrics.false_tracks = metrics.false_tracks + 1;
    end
end

%% ── 输出摘要 ─────────────────────────────────────────────────────
if ~quiet
    fprintf('\n===== 跟踪性能评估 =====\n');
    fprintf('  已确认轨迹数:  %d\n',     metrics.confirmed_tracks);
    fprintf('  虚假轨迹数:    %d\n',     metrics.false_tracks);
    fprintf('  平均确认延迟:  %.1f 帧\n', metrics.confirm_delay);
    if ~isnan(metrics.pos_rmse)
        fprintf('  位置 RMSE:     %.2f m\n', metrics.pos_rmse);
    end
    fprintf('=========================\n\n');
end
end
