%% VERIFY_EVAL_METRICS  对比原版与修正版 eval_metrics 的 RMSE 统计
clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'Evaluation'));

load(fullfile(root, 'ablation_results.mat'), 'results');
load(fullfile(root, 'Scenario', 'Step1_Data.mat'), 'truth', 'time');

tracks = results(1).tracks;
cfg    = config_tracker();

N_truth = length(truth);
N_frame = length(time);

%% 已确认轨迹索引（与 eval_metrics 一致）
conf_mask = false(1, length(tracks));
for i = 1 : length(tracks)
    if strcmp(tracks(i).status, 'confirmed') || tracks(i).total_hits >= cfg.mn_M
        conf_mask(i) = true;
    end
end
conf_idx = find(conf_mask);

%% ── 原版：真值外推 ───────────────────────────────────────────────────
truth_xy_orig = cell(N_truth, 1);
extrap_frames_per_target = zeros(N_truth, 1);
for j = 1 : N_truth
    raw = truth(j).pos3D;
    px  = raw(:, 1)';
    py  = raw(:, 2)';
    nan_mask = isnan(px) | isnan(py);
    extrap_frames_per_target(j) = sum(nan_mask);
    if any(nan_mask)
        idx = find(~nan_mask);
        px  = interp1(idx, px(idx), 1:length(px), 'linear', 'extrap');
        py  = interp1(idx, py(idx), 1:length(py), 'linear', 'extrap');
    end
    truth_xy_orig{j} = [px; py];
end
total_extrap_frames = sum(extrap_frames_per_target);

orig = collect_rmse_samples(tracks, truth, truth_xy_orig, conf_idx, N_frame, N_truth, 'original');

%% ── 修正版：仅活跃真值帧 ─────────────────────────────────────────────
fixed = collect_rmse_samples(tracks, truth, [], conf_idx, N_frame, N_truth, 'fixed');

%% ── 输出 ─────────────────────────────────────────────────────────────
fprintf('\n========== 1. 原版 eval_metrics ==========\n');
fprintf('  RMSE样本数:           %d\n', orig.n_samples);
fprintf('  有效truth数量(每帧均值): %.2f (每帧参与配对 truth 数之和 / 帧数)\n', orig.mean_truth_per_frame);
fprintf('  出现NaN后被外推的帧数: %d (跨全部目标累计)\n', total_extrap_frames);
fprintf('  NaN Truth参与评估数:   %d\n', orig.nan_truth_pairs);
fprintf('  RMSE:                  %.4f m\n', orig.rmse);

fprintf('\n========== 2. 修正版 eval_metrics ==========\n');
fprintf('  RMSE样本数:           %d\n', fixed.n_samples);
fprintf('  有效truth数量(每帧均值): %.2f\n', fixed.mean_truth_per_frame);
fprintf('  被跳过的NaN帧数:       %d (跨全部目标累计)\n', fixed.skipped_nan_frames);
fprintf('  RMSE:                  %.4f m\n', fixed.rmse);

fprintf('\n========== 3. 对比表 ==========\n');
fprintf('| 项目             | 原版      | 修正版    |\n');
fprintf('|------------------|-----------|----------|\n');
fprintf('| RMSE样本数        | %9d | %9d |\n', orig.n_samples, fixed.n_samples);
fprintf('| Truth配对数       | %9d | %9d |\n', orig.n_samples, fixed.n_samples);
fprintf('| NaN Truth参与评估数 | %9d | %9d |\n', orig.nan_truth_pairs, fixed.nan_truth_pairs);
fprintf('| RMSE (m)         | %9.2f | %9.2f |\n', orig.rmse, fixed.rmse);

%% ── Top-5 大误差样本（原版）──────────────────────────────────────────
dist_all = arrayfun(@(s) s.distance, orig.samples);
[~, ord] = sort(dist_all, 'descend');
top5 = orig.samples(ord(1:min(5, numel(ord))));

fprintf('\n========== 4. 原版 RMSE 最大的 5 个样本 ==========\n');
fprintf('%-6s %-9s %-9s %-10s %-20s\n', 'frame', 'truth_id', 'track_id', 'distance', '来源');
inactive_count = 0;
for k = 1 : length(top5)
    s = top5(k);
    reason = classify_truth_frame(truth(s.truth_id), s.frame);
    fprintf('%-6d %-9d %-9d %-10.2f %s\n', s.frame, s.truth_id, s.track_id, s.distance, reason);
    if ~strcmp(reason, '正常活跃目标')
        inactive_count = inactive_count + 1;
    end
end

pct = 100 * inactive_count / length(top5);
fprintf('\n非活跃目标占比: %d/%d (%.0f%%)\n', inactive_count, length(top5), pct);
if pct > 80
    fprintf('结论: 17.02m主要由真值外推造成。\n');
else
    fprintf('结论: 大误差样本未超过80%%来自非活跃目标，需结合其他因素分析。\n');
end

%% ══════════════════════════════════════════════════════════════════════
function out = collect_rmse_samples(tracks, truth, truth_xy_orig, conf_idx, N_frame, N_truth, mode)
samples = struct('frame', {}, 'truth_id', {}, 'track_id', {}, 'distance', {}, 'truth_active', {});
nan_truth_pairs = 0;
skipped_nan_frames = 0;
truth_per_frame_sum = 0;
n_frames_with_pairs = 0;

for k = 1 : N_frame
    est_pos = zeros(2, length(conf_idx));
    valid_est = false(1, length(conf_idx));
    for ci = 1 : length(conf_idx)
        i = conf_idx(ci);
        t_hist = tracks(i).time_hist;
        [~, loc] = min(abs(t_hist - k));
        if ~isempty(loc) && abs(t_hist(loc) - k) <= 1
            est_pos(:, ci) = tracks(i).state_hist(1:2, loc);
            valid_est(ci) = true;
        end
    end

    if strcmp(mode, 'original')
        gt_pos = zeros(2, N_truth);
        active_mask = false(1, N_truth);
        for j = 1 : N_truth
            if k <= size(truth_xy_orig{j}, 2)
                gt_pos(:, j) = truth_xy_orig{j}(:, k);
            end
            raw = truth(j).pos3D;
            active_mask(j) = k <= size(raw, 1) && ~any(isnan(raw(k, 1:2)));
        end
        n_gt = N_truth;
    else
        gt_pos = zeros(2, 0);
        active_mask = [];
        active_ids = [];
        for j = 1 : N_truth
            raw = truth(j).pos3D;
            if k <= size(raw, 1) && ~any(isnan(raw(k, 1:2)))
                gt_pos(:, end+1) = raw(k, 1:2)'; %#ok<AGROW>
                active_ids(end+1) = j; %#ok<AGROW>
            else
                skipped_nan_frames = skipped_nan_frames + 1;
            end
        end
        n_gt = size(gt_pos, 2);
        if n_gt == 0
            continue
        end
    end

    est_used = false(1, length(conf_idx));
    n_pairs_this_frame = 0;

    for j = 1 : n_gt
        if strcmp(mode, 'fixed')
            tid = active_ids(j);
            gtp = gt_pos(:, j);
            is_active = true;
        else
            tid = j;
            gtp = gt_pos(:, j);
            is_active = active_mask(j);
        end

        best_ci = -1;
        best_d  = inf;
        for ci = 1 : length(conf_idx)
            if ~valid_est(ci) || est_used(ci), continue; end
            d = norm(est_pos(:, ci) - gtp);
            if d < best_d
                best_d  = d;
                best_ci = ci;
            end
        end

        if best_ci > 0 && best_d < 100
            s.frame      = k;
            s.truth_id   = tid;
            s.track_id   = conf_idx(best_ci);
            s.distance   = best_d;
            s.truth_active = is_active;
            samples(end+1) = s; %#ok<AGROW>
            est_used(best_ci) = true;
            n_pairs_this_frame = n_pairs_this_frame + 1;
            if strcmp(mode, 'original') && ~is_active
                nan_truth_pairs = nan_truth_pairs + 1;
            end
        end
    end

    if n_pairs_this_frame > 0
        truth_per_frame_sum = truth_per_frame_sum + n_pairs_this_frame;
        n_frames_with_pairs = n_frames_with_pairs + 1;
    end
end

distances = [samples.distance];
if isempty(distances)
    rmse = NaN;
else
    rmse = sqrt(mean(distances.^2));
end

out.n_samples            = numel(distances);
out.nan_truth_pairs      = nan_truth_pairs;
out.skipped_nan_frames   = skipped_nan_frames;
out.mean_truth_per_frame = truth_per_frame_sum / max(n_frames_with_pairs, 1);
out.rmse                 = rmse;
out.samples              = samples;
end

function reason = classify_truth_frame(truth_j, frame_k)
raw = truth_j.pos3D;
if frame_k > size(raw, 1)
    reason = '目标已经消失(帧索引超范围)';
    return
end
if any(isnan(raw(frame_k, 1:2)))
    % 判断尚未出生还是已经消失
    valid_idx = find(~isnan(raw(:, 1)));
    if isempty(valid_idx)
        reason = '目标尚未出生(无有效真值)';
        return
    end
    first_v = valid_idx(1);
    last_v  = valid_idx(end);
    if frame_k < first_v
        reason = '目标尚未出生';
    elseif frame_k > last_v
        reason = '目标已经消失';
    else
        reason = '正常活跃目标(NaN间隙)';
    end
else
    reason = '正常活跃目标';
end
end
