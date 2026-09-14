%% PLOT_STEP2_DETECTION_ANALYSIS
%  Step2 点迹空间分布 + 逐帧误差/虚警统计（雷达 / 声学 / 协同）
%
%  输出目录：Scenario/figures/step2_analysis/
%
%  用法：在 Scenario/ 下运行本脚本（需 Step1_Data.mat + Step2_HeteroDetections.mat）

function plot_step2_detection_analysis()
clear; clc; close all;

this_dir    = fileparts(mfilename('fullpath'));
tracker_dir = fullfile(this_dir, '..');
out_dir     = fullfile(this_dir, 'figures', 'step2_analysis');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ── 加载数据 ─────────────────────────────────────────────────────────
s1 = fullfile(this_dir, 'Step1_Data.mat');
s2 = fullfile(this_dir, 'Step2_HeteroDetections.mat');
if ~exist(s1, 'file'), error('缺少 Step1_Data.mat'); end
if ~exist(s2, 'file'), error('缺少 Step2_HeteroDetections.mat'); end

S1 = load(s1);
S2 = load(s2, 'allDetections', 'time');

truth         = S1.truth;
time          = S2.time;
allDetections = S2.allDetections;
buildings     = S1.buildings;
acousticPos   = S1.acousticPos;
radarPos      = S1.radarPos;
if isfield(S1, 'Map'), Map = S1.Map; else, Map.size = 1000; end
if isfield(S1, 'UAV'), UAV = S1.UAV; else, UAV.T = time(end) + time(2); end

N     = length(time);
N_uav = length(truth);
colors_uav = lines(max(N_uav, 1));

% 传感器专属关联门限（米）：点到最近活跃真值 < 门限 → 判为有效探测
gate_radar    = 5.0;
gate_acoustic = 30.0;

fprintf('===== Step2 探测分析绘图 =====\n');
fprintf('  帧数: %d | 真值 UAV: %d\n', N, N_uav);
fprintf('  关联门限: 雷达 %.1f m, 声学 %.1f m\n', gate_radar, gate_acoustic);

%% ── 解析全部点迹 + 逐帧误差统计 ───────────────────────────────────────
modes = {'Radar', 'Acoustic', 'Fused'};
stats = struct();
for mi = 1:length(modes)
    stats.(modes{mi}).n_det_frame   = zeros(1, N);
    stats.(modes{mi}).n_tp_frame    = zeros(1, N);
    stats.(modes{mi}).n_fa_frame    = zeros(1, N);
    stats.(modes{mi}).rmse_frame    = nan(1, N);
    stats.(modes{mi}).mean_err_frame = nan(1, N);
    stats.(modes{mi}).all_err       = [];
    stats.(modes{mi}).all_is_tp     = logical([]);
end

rd_all = zeros(0, 2);
ac_all = zeros(0, 2);

for t = 1:N
    truth_xy = active_truth_xy(truth, t);

    [rd_xy, ac_xy, rd_err, rd_tp, ac_err, ac_tp] = ...
        analyze_frame(allDetections{t}, truth_xy, gate_radar, gate_acoustic);

    if ~isempty(rd_xy), rd_all = [rd_all; rd_xy]; end %#ok<AGROW>
    if ~isempty(ac_xy), ac_all = [ac_all; ac_xy]; end %#ok<AGROW>

    % ── 雷达 ──
    stats.Radar.n_det_frame(t) = numel(rd_tp);
    stats.Radar.n_tp_frame(t)  = sum(rd_tp);
    stats.Radar.n_fa_frame(t)  = sum(~rd_tp);
    stats.Radar = append_frame_stats(stats.Radar, rd_err, rd_tp, t);

    % ── 声学 ──
    stats.Acoustic.n_det_frame(t) = numel(ac_tp);
    stats.Acoustic.n_tp_frame(t)  = sum(ac_tp);
    stats.Acoustic.n_fa_frame(t)  = sum(~ac_tp);
    stats.Acoustic = append_frame_stats(stats.Acoustic, ac_err, ac_tp, t);

    % ── 协同（并集，各点用自身门限）──
    fused_err = [rd_err(:)', ac_err(:)'];
    fused_tp  = [rd_tp(:)', ac_tp(:)'];
    stats.Fused.n_det_frame(t) = numel(fused_tp);
    stats.Fused.n_tp_frame(t)  = sum(fused_tp);
    stats.Fused.n_fa_frame(t)  = sum(~fused_tp);
    stats.Fused = append_frame_stats(stats.Fused, fused_err, fused_tp, t);
end

% 汇总指标
summary = struct();
for mi = 1:length(modes)
    m = modes{mi};
    s = stats.(m);
    n_tot = sum(s.n_det_frame);
    n_tp  = sum(s.n_tp_frame);
    n_fa  = sum(s.n_fa_frame);
    summary.(m).n_total   = n_tot;
    summary.(m).n_tp      = n_tp;
    summary.(m).n_fa      = n_fa;
    summary.(m).fa_rate   = n_fa / max(n_tot, 1);
    summary.(m).tp_rate   = n_tp / max(n_tot, 1);
    tp_err = s.all_err(s.all_is_tp);
    if ~isempty(tp_err)
        summary.(m).rmse_all  = sqrt(mean(tp_err.^2));
        summary.(m).mean_err  = mean(tp_err);
        summary.(m).median_err = median(tp_err);
        summary.(m).p90_err   = prctile(tp_err, 90);
    else
        summary.(m).rmse_all = NaN;
        summary.(m).mean_err = NaN;
        summary.(m).median_err = NaN;
        summary.(m).p90_err = NaN;
    end
    summary.(m).fa_per_frame = mean(s.n_fa_frame);
    summary.(m).det_per_frame = mean(s.n_det_frame);
end

% 按目标计覆盖：每帧每目标最近探测误差（协同互补性）
[target_cov, n_active_frame] = target_coverage_stats( ...
    allDetections, truth, time, gate_radar, gate_acoustic);

fprintf('\n── 汇总统计 ──\n');
fprintf('%-10s %8s %8s %8s %10s %10s %10s\n', ...
    '模式', '总点迹', 'TP', 'FA', '虚警率', 'TP-RMSE', 'FA/帧');
for mi = 1:length(modes)
    m = modes{mi}; sm = summary.(m);
    fprintf('%-10s %8d %8d %8d %9.1f%% %9.2f m %9.2f\n', ...
        m, sm.n_total, sm.n_tp, sm.n_fa, 100*sm.fa_rate, sm.rmse_all, sm.fa_per_frame);
end

%% =========================================================================
%% Fig 1 — 城市地图 + 全部点迹 + 真值轨迹
%% =========================================================================
fig1 = figure('Color', 'w', 'Position', [40 40 1200 1000], 'Visible', 'off');
plot_buildings_2d(buildings, Map);
hold on; grid on; axis equal; box on;
title('Step2: All Detections on Urban Map', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

for j = 1:N_uav
    p = truth(j).pos3D;
    v = ~isnan(p(:,1));
    if sum(v) > 1
        plot(p(v,1), p(v,2), '--', 'Color', colors_uav(j,:), 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end
h_gt = plot(nan, nan, 'k--', 'LineWidth', 1.4, 'DisplayName', 'Ground Truth');

if ~isempty(ac_all)
    scatter(ac_all(:,1), ac_all(:,2), 6, [0.15 0.45 0.85], 'filled', ...
        'MarkerFaceAlpha', 0.12, 'DisplayName', ...
        sprintf('Acoustic (%d)', size(ac_all,1)));
end
if ~isempty(rd_all)
    scatter(rd_all(:,1), rd_all(:,2), 10, [0.9 0.15 0.15], 'x', ...
        'MarkerEdgeAlpha', 0.25, 'DisplayName', ...
        sprintf('Radar (%d)', size(rd_all,1)));
end
if ~isempty(acousticPos)
    scatter(acousticPos(:,1), acousticPos(:,2), 50, 'b', 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8, 'DisplayName', 'Acoustic Sensor');
end
if ~isempty(radarPos)
    scatter(radarPos(:,1), radarPos(:,2), 120, 'r', '^', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Radar');
end
legend('Location', 'northeast', 'FontSize', 9);
xlim([-0.05*Map.size, 1.05*Map.size]);
ylim([-0.05*Map.size, 1.05*Map.size]);
saveas(fig1, fullfile(out_dir, 'fig01_map_all_detections.png'));
close(fig1);
fprintf('  已保存 fig01_map_all_detections.png\n');

%% Fig 1b — 分传感器子图
fig1b = figure('Color', 'w', 'Position', [60 60 1300 560], 'Visible', 'off');
titles = {'Radar Only', 'Acoustic Only', 'Fused (Overlay)'};
datasets = {rd_all, ac_all, []};
for sp = 1:3
    subplot(1, 3, sp);
    plot_buildings_2d(buildings, Map);
    hold on; grid on; axis equal; box on;
    title(titles{sp}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('X (m)'); ylabel('Y (m)');
    for j = 1:N_uav
        p = truth(j).pos3D; v = ~isnan(p(:,1));
        if sum(v) > 1
            plot(p(v,1), p(v,2), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    if sp == 1 && ~isempty(rd_all)
        scatter(rd_all(:,1), rd_all(:,2), 8, [0.9 0.2 0.2], 'x', 'MarkerEdgeAlpha', 0.3);
    elseif sp == 2 && ~isempty(ac_all)
        scatter(ac_all(:,1), ac_all(:,2), 6, [0.2 0.45 0.9], 'filled', 'MarkerFaceAlpha', 0.15);
    elseif sp == 3
        if ~isempty(ac_all)
            scatter(ac_all(:,1), ac_all(:,2), 5, [0.2 0.45 0.9], 'filled', ...
                'MarkerFaceAlpha', 0.1, 'DisplayName', 'Acoustic');
        end
        if ~isempty(rd_all)
            scatter(rd_all(:,1), rd_all(:,2), 8, [0.9 0.2 0.2], 'x', ...
                'MarkerEdgeAlpha', 0.3, 'DisplayName', 'Radar');
        end
        legend('Location', 'best', 'FontSize', 8);
    end
    if ~isempty(acousticPos)
        scatter(acousticPos(:,1), acousticPos(:,2), 28, 'b', 'filled', 'HandleVisibility', 'off');
    end
    if ~isempty(radarPos)
        scatter(radarPos(:,1), radarPos(:,2), 70, 'r', '^', 'filled', 'HandleVisibility', 'off');
    end
    xlim([-0.05*Map.size, 1.05*Map.size]);
    ylim([-0.05*Map.size, 1.05*Map.size]);
end
sgtitle('Step2 Detections by Sensor Mode', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig1b, fullfile(out_dir, 'fig02_map_by_sensor.png'));
close(fig1b);
fprintf('  已保存 fig02_map_by_sensor.png\n');

%% =========================================================================
%% Fig 2 — 逐帧 RMSE（仅 TP）与虚警数量
%% =========================================================================
col_rd = [0.85 0.25 0.25];
col_ac = [0.25 0.45 0.85];
col_fu = [0.15 0.65 0.35];

fig2 = figure('Color', 'w', 'Position', [80 80 1200 820], 'Visible', 'off');

subplot(2, 1, 1);
hold on; grid on; box on;
plot(time, stats.Radar.rmse_frame, '-', 'Color', col_rd, 'LineWidth', 1.2, 'DisplayName', 'Radar');
plot(time, stats.Acoustic.rmse_frame, '-', 'Color', col_ac, 'LineWidth', 1.2, 'DisplayName', 'Acoustic');
plot(time, stats.Fused.rmse_frame, '-', 'Color', col_fu, 'LineWidth', 1.4, 'DisplayName', 'Fused');
title('Per-Frame RMSE (True Positives Only)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('RMSE (m)');
legend('Location', 'northeast');
xlim([0, time(end)]);

subplot(2, 1, 2);
hold on; grid on; box on;
plot(time, stats.Radar.n_fa_frame, '-', 'Color', col_rd, 'LineWidth', 1.2, 'DisplayName', 'Radar FA');
plot(time, stats.Acoustic.n_fa_frame, '-', 'Color', col_ac, 'LineWidth', 1.2, 'DisplayName', 'Acoustic FA');
plot(time, stats.Fused.n_fa_frame, '-', 'Color', col_fu, 'LineWidth', 1.4, 'DisplayName', 'Fused FA');
plot(time, stats.Radar.n_tp_frame, '--', 'Color', col_rd*0.7+0.3, 'LineWidth', 0.9, 'DisplayName', 'Radar TP');
plot(time, stats.Acoustic.n_tp_frame, '--', 'Color', col_ac*0.7+0.3, 'LineWidth', 0.9, 'DisplayName', 'Acoustic TP');
title('Per-Frame True Positive / False Alarm Count', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Count');
legend('Location', 'northeast', 'FontSize', 8);
xlim([0, time(end)]);

sgtitle('Step2: Frame-wise Error and False Alarm', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig2, fullfile(out_dir, 'fig03_per_frame_rmse_fa.png'));
close(fig2);
fprintf('  已保存 fig03_per_frame_rmse_fa.png\n');

%% Fig 3 — 逐帧平均最近误差（含 FA 的大误差）
fig3 = figure('Color', 'w', 'Position', [90 90 1200 500], 'Visible', 'off');
hold on; grid on; box on;
mean_nearest = nan(3, N);
for t = 1:N
    truth_xy = active_truth_xy(truth, t);
    [~, ~, rd_e, ~, ac_e, ~] = analyze_frame(allDetections{t}, truth_xy, gate_radar, gate_acoustic);
    mean_nearest(1, t) = mean_or_nan(rd_e);
    mean_nearest(2, t) = mean_or_nan(ac_e);
    mean_nearest(3, t) = mean_or_nan([rd_e, ac_e]);
end
plot(time, mean_nearest(1,:), '-', 'Color', col_rd, 'LineWidth', 1.2, 'DisplayName', 'Radar');
plot(time, mean_nearest(2,:), '-', 'Color', col_ac, 'LineWidth', 1.2, 'DisplayName', 'Acoustic');
plot(time, mean_nearest(3,:), '-', 'Color', col_fu, 'LineWidth', 1.4, 'DisplayName', 'Fused');
yline(gate_radar, ':', 'Color', col_rd, 'LineWidth', 1, 'DisplayName', sprintf('Radar gate %.0f m', gate_radar));
yline(gate_acoustic, ':', 'Color', col_ac, 'LineWidth', 1, 'DisplayName', sprintf('Acoustic gate %.0f m', gate_acoustic));
title('Per-Frame Mean Distance to Nearest Truth (All Detections)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Mean nearest distance (m)');
legend('Location', 'northeast', 'FontSize', 8);
xlim([0, time(end)]);
saveas(fig3, fullfile(out_dir, 'fig04_per_frame_mean_nearest.png'));
close(fig3);
fprintf('  已保存 fig04_per_frame_mean_nearest.png\n');

%% =========================================================================
%% Fig 4 — 统计汇总：虚警率、误差分布、CDF
%% =========================================================================
fig4 = figure('Color', 'w', 'Position', [100 100 1300 900], 'Visible', 'off');

% 4a 虚警率 / 检测率柱状图
subplot(2, 2, 1);
fa_rates = [summary.Radar.fa_rate, summary.Acoustic.fa_rate, summary.Fused.fa_rate] * 100;
tp_rates = [summary.Radar.tp_rate, summary.Acoustic.tp_rate, summary.Fused.tp_rate] * 100;
b = bar([fa_rates; tp_rates]', 'grouped');
b(1).FaceColor = [0.9 0.35 0.35];
b(2).FaceColor = [0.35 0.75 0.45];
set(gca, 'XTickLabel', modes);
ylabel('Rate (%)');
title('False Alarm vs True Positive Rate', 'FontWeight', 'bold');
legend({'False Alarm', 'True Positive'}, 'Location', 'northwest');
grid on;

% 4b TP 误差箱线图
subplot(2, 2, 2);
tp_data = {stats.Radar.all_err(stats.Radar.all_is_tp), ...
           stats.Acoustic.all_err(stats.Acoustic.all_is_tp), ...
           stats.Fused.all_err(stats.Fused.all_is_tp)};
boxplot([tp_data{1}(:); tp_data{2}(:); tp_data{3}(:)], ...
    [ones(numel(tp_data{1}),1); 2*ones(numel(tp_data{2}),1); 3*ones(numel(tp_data{3}),1)], ...
    'Labels', modes, 'Whisker', 1.5);
ylabel('Error to nearest truth (m)');
title('TP Error Distribution (Boxplot)', 'FontWeight', 'bold');
grid on;

% 4c 误差直方图
subplot(2, 2, 3);
hold on;
histogram(tp_data{1}, 40, 'FaceColor', col_rd, 'FaceAlpha', 0.45, 'EdgeColor', 'none', 'DisplayName', 'Radar');
histogram(tp_data{2}, 40, 'FaceColor', col_ac, 'FaceAlpha', 0.45, 'EdgeColor', 'none', 'DisplayName', 'Acoustic');
xline(gate_radar, '--', 'Color', col_rd, 'LineWidth', 1.2);
xline(gate_acoustic, '--', 'Color', col_ac, 'LineWidth', 1.2);
xlabel('Error (m)'); ylabel('Count');
title('TP Error Histogram', 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on;

% 4d CDF
subplot(2, 2, 4);
hold on;
for mi = 1:3
    e = tp_data{mi};
    if ~isempty(e)
        [fx, x] = ecdf(e);
        cols = {col_rd, col_ac, col_fu};
        plot(x, fx, 'LineWidth', 1.8, 'Color', cols{mi}, 'DisplayName', modes{mi});
    end
end
xlabel('Error (m)'); ylabel('CDF');
title('TP Error CDF', 'FontWeight', 'bold');
legend('Location', 'southeast');
grid on;

sgtitle('Step2: Detection Error Statistics', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig4, fullfile(out_dir, 'fig05_error_statistics.png'));
close(fig4);
fprintf('  已保存 fig05_error_statistics.png\n');

%% Fig 5 — 按目标覆盖：每帧最近探测误差（协同互补）
fig5 = figure('Color', 'w', 'Position', [110 110 1200 520], 'Visible', 'off');
hold on; grid on; box on;
plot(time, target_cov.radar, '-', 'Color', col_rd, 'LineWidth', 1.2, 'DisplayName', 'Radar best/target');
plot(time, target_cov.acoustic, '-', 'Color', col_ac, 'LineWidth', 1.2, 'DisplayName', 'Acoustic best/target');
plot(time, target_cov.fused_best, '-', 'Color', col_fu, 'LineWidth', 1.6, 'DisplayName', 'Fused min(radar,acoustic)');
title('Per-Frame Best Detection Error per Active Target', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Mean best error (m)');
legend('Location', 'northeast');
xlim([0, time(end)]);
subtitle(sprintf('Active frames: %d | Fused improves over radar-only in %.1f%% of active frames', ...
    n_active_frame, 100 * target_cov.fused_better_frac));
saveas(fig5, fullfile(out_dir, 'fig06_target_best_error.png'));
close(fig5);
fprintf('  已保存 fig06_target_best_error.png\n');

%% Fig 6 — 虚警率对比 + 每帧探测密度
fig6 = figure('Color', 'w', 'Position', [120 120 1200 500], 'Visible', 'off');
subplot(1, 2, 1);
fa_mean_frame = mean([stats.Radar.n_fa_frame; stats.Acoustic.n_fa_frame; stats.Fused.n_fa_frame], 2);
bar(fa_mean_frame, 'FaceColor', 'flat', 'CData', [col_rd; col_ac; col_fu]);
set(gca, 'XTickLabel', modes);
ylabel('Mean FA count per frame');
title('Average False Alarms per Frame', 'FontWeight', 'bold');
grid on;
for i = 1:3
    text(i, fa_mean_frame(i) + 0.3, sprintf('%.1f%% FA', fa_rates(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

subplot(1, 2, 2);
bar(mean([stats.Radar.n_det_frame; stats.Acoustic.n_det_frame; stats.Fused.n_det_frame], 2), ...
    'FaceColor', 'flat', 'CData', [col_rd; col_ac; col_fu]);
set(gca, 'XTickLabel', modes);
ylabel('Mean detections per frame');
title('Average Detection Load per Frame', 'FontWeight', 'bold');
grid on;

sgtitle('False Alarm Load Comparison', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig6, fullfile(out_dir, 'fig07_fa_load_comparison.png'));
close(fig6);
fprintf('  已保存 fig07_fa_load_comparison.png\n');

%% 保存数值结果
results = struct('summary', summary, 'stats', stats, 'target_cov', target_cov, ...
    'gate_radar', gate_radar, 'gate_acoustic', gate_acoustic, 'modes', {modes});
save(fullfile(out_dir, 'step2_analysis_results.mat'), 'results', 'summary', 'target_cov');
fprintf('\n===== 分析完成 | 输出: %s =====\n', out_dir);

end

%% ══════════════════════════════════════════════════════════════════════
function s = append_frame_stats(s, errs, is_tp, t)
if isempty(errs) || isempty(is_tp) || numel(errs) ~= numel(is_tp)
    return
end
errs  = errs(:)';
is_tp = logical(is_tp(:)');
tp_e  = errs(is_tp);
fa_e  = errs(~is_tp);
s.all_err   = [s.all_err, errs]; %#ok<AGROW>
s.all_is_tp = [s.all_is_tp, is_tp]; %#ok<AGROW>
if ~isempty(tp_e)
    s.rmse_frame(t)     = sqrt(mean(tp_e.^2));
    s.mean_err_frame(t) = mean(tp_e);
end
end

function truth_xy = active_truth_xy(truth, t)
truth_xy = zeros(2, 0);
for j = 1:length(truth)
    p = truth(j).pos3D;
    if t <= size(p, 1) && ~any(isnan(p(t, 1:2)))
        truth_xy(:, end+1) = p(t, 1:2)'; %#ok<AGROW>
    end
end
end

function [rd_xy, ac_xy, rd_err, rd_tp, ac_err, ac_tp] = analyze_frame( ...
    frame, truth_xy, gate_radar, gate_acoustic)
rd_xy = zeros(0, 2); ac_xy = zeros(0, 2);
rd_err = []; rd_tp = [];
ac_err = []; ac_tp = [];
if isempty(frame) || isempty(truth_xy), return; end

if iscell(frame), det_list = frame; else, det_list = num2cell(frame); end
n_truth = size(truth_xy, 2);

for d = 1:length(det_list)
    det = det_list{d};
    if isobject(det)
        meas = det.Measurement(:);
        sid  = det.SensorIndex;
    elseif isstruct(det)
        meas = det.Measurement(:);
        sid  = det.SensorIndex;
    else
        continue
    end
    if numel(meas) < 2, continue; end
    z = meas(1:2);

    dists = sqrt(sum((truth_xy - z).^2, 1));
    min_d = min(dists);

    if sid > 100
        gate = gate_radar;
        rd_xy = [rd_xy; z']; %#ok<AGROW>
        rd_err(end+1) = min_d; %#ok<AGROW>
        rd_tp(end+1)  = min_d <= gate; %#ok<AGROW>
    else
        gate = gate_acoustic;
        ac_xy = [ac_xy; z']; %#ok<AGROW>
        ac_err(end+1) = min_d; %#ok<AGROW>
        ac_tp(end+1)  = min_d <= gate; %#ok<AGROW>
    end
end
end

function [cov, n_active] = target_coverage_stats(allDetections, truth, time, gate_radar, gate_acoustic)
N = length(time);
cov.radar = nan(1, N);
cov.acoustic = nan(1, N);
cov.fused_best = nan(1, N);
n_better = 0; n_active = 0;

for t = 1:N
    truth_xy = active_truth_xy(truth, t);
    if isempty(truth_xy), continue; end
    n_active = n_active + 1;
    n_tgt = size(truth_xy, 2);

    [rd_xy, ac_xy, ~, ~, ~, ~] = analyze_frame(allDetections{t}, truth_xy, gate_radar, gate_acoustic);

    rd_best = nan(1, n_tgt);
    ac_best = nan(1, n_tgt);
    for j = 1:n_tgt
        tx = truth_xy(:, j);
        if ~isempty(rd_xy)
            rd_best(j) = min(sqrt(sum((rd_xy' - tx).^2, 1)));
        end
        if ~isempty(ac_xy)
            ac_best(j) = min(sqrt(sum((ac_xy' - tx).^2, 1)));
        end
    end

    cov.radar(t)    = mean(rd_best, 'omitnan');
    cov.acoustic(t) = mean(ac_best, 'omitnan');

    fused = nan(1, n_tgt);
    for j = 1:n_tgt
        candidates = [rd_best(j), ac_best(j)];
        fused(j) = min(candidates, [], 'omitnan');
    end
    cov.fused_best(t) = mean(fused, 'omitnan');

    if ~isnan(cov.radar(t)) && ~isnan(cov.fused_best(t)) && cov.fused_best(t) < cov.radar(t) - 0.5
        n_better = n_better + 1;
    end
end
cov.fused_better_frac = n_better / max(n_active, 1);
end

function v = mean_or_nan(x)
if isempty(x), v = NaN; else, v = mean(x); end
end

function plot_buildings_2d(buildings, Map)
hold on;
maxH = 1;
for i = 1:length(buildings)
    b = buildings(i);
    if isfield(b, 'pos')
        bx = b.pos(1); by = b.pos(2);
        bw = b.dim(1); bd = b.dim(2); bh = b.dim(3);
    else
        bx = b.center(1) - b.dim(1)/2;
        by = b.center(2) - b.dim(2)/2;
        bw = b.dim(1); bd = b.dim(2); bh = b.dim(3);
    end
    maxH = max(maxH, bh);
    fill([bx, bx+bw, bx+bw, bx], [by, by, by+bd, by+bd], ...
        [0.55 0.55 0.6], 'FaceAlpha', 0.25 + 0.4*(bh/maxH), ...
        'EdgeColor', [0.35 0.35 0.4], 'LineWidth', 0.4, 'HandleVisibility', 'off');
end
end
