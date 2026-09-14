%% =========================================================================
%  plot_full_pipeline.m
%  全流程可视化：Step1 真值/传感器 → Step2 异构探测 → Step3 跟踪输出
%
%  输入（优先 Scenario/ 目录）：
%    Step1_Data.mat, Step2_HeteroDetections.mat, Step3_IMM_JPDA_EKF_Result.mat
%
%  输出目录：Scenario/figures/
%    fig01~04  若不存在则跳过（由 run_scenario 生成）
%    fig05     Step2 俯视图：真值 + 雷达/声学点迹 + 传感器
%    fig06     Step2 每帧探测数量时序
%    fig07     Step2 典型帧点迹快照（4 子图）
%    fig08     Step3 俯视图：真值 vs 估计轨迹（已确认）
%    fig09     Step3 轨迹维护：确认/暂定数量随时间
%    fig10     Step3 RMSE 与 IMM 模型概率
%    fig11     Step3 轨迹生命周期（积分式维护：诞生→确认→删除）
%    fig12     三阶段对照（真值 / 观测 / 跟踪）一张图看清全流程
%% =========================================================================
function plot_full_pipeline()
clear; clc; close all;

this_dir   = fileparts(mfilename('fullpath'));
tracker_dir = fullfile(this_dir, '..');
v2_dir     = fullfile(tracker_dir, '..', 'v2.0');
out_dir    = fullfile(this_dir, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

addpath(this_dir);
addpath(tracker_dir);
addpath(fullfile(tracker_dir, 'Evaluation'));
addpath(v2_dir);

cfg = config_tracker();

%% ── 加载 Step1 ───────────────────────────────────────────────────────
s1 = fullfile(this_dir, 'Step1_Data.mat');
if ~exist(s1, 'file')
    error('缺少 %s，请先运行 run_scenario.m', s1);
end
S1 = load(s1);
truth      = S1.truth;
time       = S1.time;
buildings  = S1.buildings;
acousticPos = S1.acousticPos;
radarPos   = S1.radarPos;
if isfield(S1, 'Map'), Map = S1.Map; else, Map.size = 1000; end
if isfield(S1, 'UAV'), UAV = S1.UAV; else, UAV.T = time(end)+time(2); end
N = length(time);
colors_uav = lines(max(1, length(truth)));

fprintf('===== 全流程可视化 =====\n');
fprintf('  Step1: %s\n', s1);

%% ── 加载 Step2 ───────────────────────────────────────────────────────
s2 = fullfile(this_dir, 'Step2_HeteroDetections.mat');
if ~exist(s2, 'file')
    s2v2 = fullfile(v2_dir, 'Step2_HeteroDetections.mat');
    if exist(s2v2, 'file')
        s2 = s2v2;
    else
        error('缺少 Step2_HeteroDetections.mat，请先运行 step2_detection_simulation.m');
    end
end
S2 = load(s2, 'allDetections', 'time');
allDetections = S2.allDetections;
if length(S2.time) ~= N, time = S2.time; N = length(time); end
fprintf('  Step2: %s\n', s2);

%% ── 加载 Step3（可选）────────────────────────────────────────────────
tracks = []; metrics = []; frame_log = [];
s3 = fullfile(tracker_dir, 'Step3_IMM_JPDA_EKF_Result.mat');
if ~exist(s3, 'file')
    s3 = fullfile(this_dir, 'Step3_IMM_JPDA_EKF_Result.mat');
end
has_step3 = exist(s3, 'file');
if has_step3
    S3 = load(s3, 'tracks', 'metrics', 'frame_log', 'cfg');
    tracks    = S3.tracks;
    metrics   = S3.metrics;
    frame_log = S3.frame_log;
    if isfield(S3, 'cfg'), cfg = S3.cfg; end
    fprintf('  Step3: %s\n', s3);
else
    fprintf('  Step3: 未找到结果，仅生成 Step2 图；运行 main_imm_jpda_ekf 后可重跑本脚本。\n');
end

%% ── 解析 Step2 点迹 ──────────────────────────────────────────────────
[rd_xy, ac_xy, n_rd, n_ac] = extract_all_detections_xy(allDetections);

%% =========================================================================
%% Fig05 — Step2 俯视图：真值 + 点迹 + 传感器
%% =========================================================================
fig5 = figure('Color','w','Position',[60 60 950 850],'Visible','off');
plot_buildings_2d(buildings, Map);
hold on; grid on; axis equal; box on;
title('Step2: Heterogeneous Detections (Top-Down)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

for i = 1:length(truth)
    p = truth(i).pos3D;
    v = ~isnan(p(:,1));
    if sum(v) > 1
        plot(p(v,1), p(v,2), '--', 'Color', colors_uav(i,:), 'LineWidth', 1.4, ...
            'HandleVisibility','off');
    end
end
h_truth = plot(nan, nan, 'k--', 'LineWidth', 1.4, 'DisplayName', 'Ground Truth');

if ~isempty(ac_xy)
    scatter(ac_xy(:,1), ac_xy(:,2), 8, [0.2 0.45 0.9], 'filled', ...
        'MarkerFaceAlpha', 0.15, 'DisplayName', 'Acoustic');
end
if ~isempty(rd_xy)
    scatter(rd_xy(:,1), rd_xy(:,2), 14, [0.9 0.2 0.2], 'x', ...
        'MarkerEdgeAlpha', 0.35, 'DisplayName', 'Radar');
end
if ~isempty(acousticPos)
    scatter(acousticPos(:,1), acousticPos(:,2), 42, 'b', 'filled', ...
        'MarkerEdgeColor','w', 'DisplayName', 'Acoustic Sensor');
end
if ~isempty(radarPos)
    scatter(radarPos(:,1), radarPos(:,2), 110, 'r', '^', 'filled', ...
        'MarkerEdgeColor','k', 'DisplayName', 'Radar');
end
legend('Location','northeast', 'FontSize', 9);
xlim([-0.05*Map.size, 1.05*Map.size]);
ylim([-0.05*Map.size, 1.05*Map.size]);
saveas(fig5, fullfile(out_dir, 'fig05_step2_detections_topdown.png'));
close(fig5);
fprintf('  已保存 fig05_step2_detections_topdown.png\n');

%% =========================================================================
%% Fig06 — 每帧探测数量
%% =========================================================================
fig6 = figure('Color','w','Position',[80 80 1000 420],'Visible','off');
hold on; grid on; box on;
ah = area(time, [n_rd(:), n_ac(:)], 'LineStyle', 'none');
ah(1).FaceColor = [0.85 0.35 0.35]; ah(1).FaceAlpha = 0.8;
ah(2).FaceColor = [0.35 0.55 0.85]; ah(2).FaceAlpha = 0.8;
plot(time, n_rd + n_ac, 'k-', 'LineWidth', 1.2);
title('Step2: Detections per Frame', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Count');
legend([ah(1), ah(2)], {'Radar', 'Acoustic'}, 'Location','northeast');
xlim([0, time(end)]);
saveas(fig6, fullfile(out_dir, 'fig06_step2_detections_per_frame.png'));
close(fig6);
fprintf('  已保存 fig06_step2_detections_per_frame.png\n');

%% =========================================================================
%% Fig07 — 典型帧快照
%% =========================================================================
snap_t = [10, 20, 30, 40];
snap_t = snap_t(snap_t <= time(end));
fig7 = figure('Color','w','Position',[90 90 1100 900],'Visible','off');
for si = 1:length(snap_t)
    subplot(2, 2, si);
    plot_buildings_2d(buildings, Map);
    hold on; grid on; axis equal; box on;
    tk = find(abs(time - snap_t(si)) < 1e-6, 1);
    if isempty(tk), tk = round(snap_t(si) / cfg.dt) + 1; end
    tk = min(max(tk, 1), N);

    for j = 1:length(truth)
        p = truth(j).pos3D;
        if tk <= size(p,1) && ~any(isnan(p(tk,:)))
            plot(p(tk,1), p(tk,2), 'o', 'Color', colors_uav(j,:), ...
                'MarkerSize', 10, 'MarkerFaceColor', colors_uav(j,:), ...
                'HandleVisibility','off');
            text(p(tk,1)+12, p(tk,2)+12, sprintf('#%d', j), ...
                'Color', colors_uav(j,:), 'FontWeight', 'bold', 'FontSize', 9);
        end
    end
    [zr, za] = frame_detections_xy(allDetections{tk});
    if ~isempty(za), scatter(za(:,1), za(:,2), 22, 'b', 'filled', 'MarkerFaceAlpha', 0.5); end
    if ~isempty(zr), scatter(zr(:,1), zr(:,2), 36, 'r', 'x', 'LineWidth', 1.2); end
    title(sprintf('t = %.1f s  (%d det)', time(tk), n_rd(tk)+n_ac(tk)), 'FontSize', 11);
    xlim([-0.05*Map.size, 1.05*Map.size]);
    ylim([-0.05*Map.size, 1.05*Map.size]);
end
sgtitle('Step2: Detection Snapshots at Key Times', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig7, fullfile(out_dir, 'fig07_step2_frame_snapshots.png'));
close(fig7);
fprintf('  已保存 fig07_step2_frame_snapshots.png\n');

if ~has_step3
    copy_step1_figures(this_dir, out_dir);
    fprintf('\n===== Step2 图已生成；补跑 main_imm_jpda_ekf 后重跑本脚本可得 Step3 图 =====\n');
    fprintf('  输出目录: %s\n', out_dir);
    return
end

%% =========================================================================
%% Fig08 — Step3 真值 vs 估计（俯视图）
%% =========================================================================
fig8 = figure('Color','w','Position',[70 70 980 860],'Visible','off');
plot_buildings_2d(buildings, Map);
hold on; grid on; axis equal; box on;
title('Step3: Ground Truth vs Track Estimates', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

for j = 1:length(truth)
    p = truth(j).pos3D;
    v = ~isnan(p(:,1));
    if sum(v) > 1
        plot(p(v,1), p(v,2), 'k--', 'LineWidth', 1.3, 'HandleVisibility','off');
    end
end
h_gt = plot(nan, nan, 'k--', 'LineWidth', 1.3, 'DisplayName', 'Ground Truth');

colors_tr = lines(max(1, length(tracks)));
shown = 0;
for i = 1:length(tracks)
    if tracks(i).total_hits < cfg.mn_M, continue; end
    sh = tracks(i).state_hist;
    c  = colors_tr(mod(i-1, size(colors_tr,1))+1, :);
    plot(sh(1,:), sh(2,:), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(sh(1,1), sh(2,1), 'o', 'Color', c, 'MarkerSize', 7, ...
        'MarkerFaceColor', c, 'HandleVisibility','off');
    text(sh(1,end)+8, sh(2,end)+8, sprintf('ID%d', tracks(i).id), ...
        'Color', c, 'FontSize', 8, 'FontWeight', 'bold');
    shown = shown + 1;
end
h_est = plot(nan, nan, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Confirmed Estimate');
legend([h_gt, h_est], 'Location', 'northeast', 'FontSize', 9);
subtitle(sprintf('%d confirmed tracks (M/N=%d/%d)', shown, cfg.mn_M, cfg.mn_N));
xlim([-0.05*Map.size, 1.05*Map.size]);
ylim([-0.05*Map.size, 1.05*Map.size]);
saveas(fig8, fullfile(out_dir, 'fig08_step3_tracking_topdown.png'));
close(fig8);
fprintf('  已保存 fig08_step3_tracking_topdown.png\n');

%% =========================================================================
%% Fig09 — 轨迹维护：帧级 confirmed / tentative 数量
%% =========================================================================
fig9 = figure('Color','w','Position',[90 90 1000 440],'Visible','off');
hold on; grid on; box on;
if ~isempty(frame_log) && isfield(frame_log, 'n_confirmed')
    plot(time, frame_log.n_confirmed, 'b-', 'LineWidth', 2, 'DisplayName', 'Confirmed');
    plot(time, frame_log.n_tentative, 'Color', [0.9 0.55 0.1], 'LineWidth', 1.5, ...
        'DisplayName', 'Tentative');
    plot(time, frame_log.n_deleted, 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'DisplayName', 'Deleted (cumulative status)');
    yline(length(truth), 'k--', sprintf('%d truth UAV', length(truth)), ...
        'LabelHorizontalAlignment', 'left');
end
title('Step3: Track Maintenance — Active Count per Frame', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Track Count');
legend('Location', 'northeast');
xlim([0, time(end)]);
saveas(fig9, fullfile(out_dir, 'fig09_step3_track_counts.png'));
close(fig9);
fprintf('  已保存 fig09_step3_track_counts.png\n');

%% =========================================================================
%% Fig10 — RMSE + IMM 模型概率
%% =========================================================================
fig10 = figure('Color','w','Position',[100 100 1100 480],'Visible','off');
subplot(1,2,1); hold on; grid on; box on;
rmse_ts = compute_rmse_timeseries(tracks, truth, time, cfg);
vk = ~isnan(rmse_ts);
if any(vk)
    plot(time(vk), rmse_ts(vk), 'b-', 'LineWidth', 1.5);
end
yline(5, 'r--', '5 m ref');
title('Step3: Per-Frame Position RMSE', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('RMSE (m)');
if isstruct(metrics) && isfield(metrics, 'pos_rmse') && ~isnan(metrics.pos_rmse)
    text(0.02, 0.95, sprintf('Mean RMSE = %.2f m', metrics.pos_rmse), ...
        'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
end

subplot(1,2,2); hold on; grid on; box on;
title('Step3: IMM Model Probabilities (sample track)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frame'); ylabel('mu_j');
plotted = false;
for i = 1:length(tracks)
    if tracks(i).total_hits >= cfg.mn_M && ~isempty(tracks(i).mu_hist)
        mu_h = tracks(i).mu_hist;
        th   = tracks(i).time_hist;
        for j = 1:size(mu_h, 1)
            plot(th, mu_h(j,:), 'LineWidth', 1.4, 'DisplayName', sprintf('Model %d', j));
        end
        subtitle(sprintf('Track ID %d', tracks(i).id));
        plotted = true;
        break
    end
end
if ~plotted
    text(0.5, 0.5, 'No confirmed track with mu history', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
end
legend('Location', 'northeast', 'FontSize', 8);
ylim([0, 1]);
saveas(fig10, fullfile(out_dir, 'fig10_step3_rmse_imm.png'));
close(fig10);
fprintf('  已保存 fig10_step3_rmse_imm.png\n');

%% =========================================================================
%% Fig11 — 积分式轨迹维护：每条航迹的诞生 / 确认 / 存活区间
%% =========================================================================
fig11 = figure('Color','w','Position',[110 110 1100 520],'Visible','off');
hold on; grid on; box on;
y_pos = 1;
yticks_lbl = {};
yticks_pos = [];

for i = 1:length(tracks)
    tr = tracks(i);
    t0 = tr.birth_frame;
    t1 = tr.time_hist(end);
    if t0 < 1 || t1 < t0, continue; end

    t_ax = [time(t0), time(t1)];
    is_conf = tr.total_hits >= cfg.mn_M;
    switch tr.status
        case 'confirmed', col = [0.15 0.55 0.25];
        case 'tentative',  col = [0.95 0.65 0.15];
        otherwise,        col = [0.55 0.55 0.55];
    end
    if is_conf && strcmp(tr.status, 'deleted')
        col = [0.35 0.65 0.45];
    end

    patch([t_ax(1) t_ax(2) t_ax(2) t_ax(1)], [y_pos-0.35 y_pos-0.35 y_pos+0.35 y_pos+0.35], ...
        col, 'FaceAlpha', 0.75, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.6);
    text(t_ax(1)+0.3, y_pos, sprintf('#%d', tr.id), 'FontSize', 8, ...
        'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

    yticks_lbl{end+1} = sprintf('T%d', tr.id); %#ok<AGROW>
    yticks_pos(end+1) = y_pos; %#ok<AGROW>
    y_pos = y_pos + 1;
end

for j = 1:length(truth)
    bt = truth(j).birthTime;
    dt = truth(j).deathTime;
    yline(bt, ':', 'Color', colors_uav(j,:), 'LineWidth', 0.9, 'HandleVisibility','off');
    yline(dt, '-.', 'Color', colors_uav(j,:), 'LineWidth', 0.9, 'HandleVisibility','off');
end

set(gca, 'YTick', yticks_pos, 'YTickLabel', yticks_lbl);
title('Step3: Incremental Track Maintenance (Birth \rightarrow Life \rightarrow Delete)', ...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Track ID');
xlim([0, time(end)]);

h_c = patch(nan, nan, [0.15 0.55 0.25], 'DisplayName', 'Confirmed');
h_t = patch(nan, nan, [0.95 0.65 0.15], 'DisplayName', 'Tentative');
h_d = patch(nan, nan, [0.55 0.55 0.55], 'DisplayName', 'Deleted');
legend([h_c, h_t, h_d], 'Location', 'northeast', 'FontSize', 9);
saveas(fig11, fullfile(out_dir, 'fig11_step3_track_lifecycle.png'));
close(fig11);
fprintf('  已保存 fig11_step3_track_lifecycle.png\n');

%% =========================================================================
%% Fig12 — 三阶段对照（全流程一图）
%% =========================================================================
fig12 = figure('Color','w','Position',[40 40 1400 480],'Visible','off');
stages = {'(a) Step1: Truth Trajectories', ...
          '(b) Step2: Sensor Observations', ...
          '(c) Step3: Filtered Track Output'};
for sp = 1:3
    subplot(1, 3, sp);
    plot_buildings_2d(buildings, Map);
    hold on; grid on; axis equal; box on;
    title(stages{sp}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('X (m)'); ylabel('Y (m)');

    switch sp
        case 1
            for i = 1:length(truth)
                p = truth(i).pos3D;
                v = ~isnan(p(:,1));
                if sum(v) > 1
                    plot(p(v,1), p(v,2), '-', 'Color', colors_uav(i,:), 'LineWidth', 1.6);
                end
            end
            if ~isempty(acousticPos)
                scatter(acousticPos(:,1), acousticPos(:,2), 20, 'b', 'filled', ...
                    'MarkerFaceAlpha', 0.6);
            end
            if ~isempty(radarPos)
                scatter(radarPos(:,1), radarPos(:,2), 60, 'r', '^', 'filled');
            end

        case 2
            for i = 1:length(truth)
                p = truth(i).pos3D;
                v = ~isnan(p(:,1));
                if sum(v) > 1
                    plot(p(v,1), p(v,2), 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
                end
            end
            if ~isempty(ac_xy)
                scatter(ac_xy(:,1), ac_xy(:,2), 6, [0.2 0.45 0.9], 'filled', ...
                    'MarkerFaceAlpha', 0.12);
            end
            if ~isempty(rd_xy)
                scatter(rd_xy(:,1), rd_xy(:,2), 10, [0.9 0.2 0.2], 'x', ...
                    'MarkerEdgeAlpha', 0.25);
            end

        case 3
            for j = 1:length(truth)
                p = truth(j).pos3D;
                v = ~isnan(p(:,1));
                if sum(v) > 1
                    plot(p(v,1), p(v,2), 'k--', 'LineWidth', 0.9, 'HandleVisibility','off');
                end
            end
            for i = 1:length(tracks)
                if tracks(i).total_hits < cfg.mn_M, continue; end
                sh = tracks(i).state_hist;
                plot(sh(1,:), sh(2,:), '-', 'LineWidth', 1.4);
            end
    end
    xlim([-0.05*Map.size, 1.05*Map.size]);
    ylim([-0.05*Map.size, 1.05*Map.size]);
end
sgtitle('Full Pipeline: Truth \rightarrow Detections \rightarrow IMM-JPDA-EKF Tracks', ...
    'FontSize', 13, 'FontWeight', 'bold');
saveas(fig12, fullfile(out_dir, 'fig12_pipeline_three_stages.png'));
close(fig12);
fprintf('  已保存 fig12_pipeline_three_stages.png\n');

copy_step1_figures(this_dir, out_dir);

fprintf('\n===== 全流程图生成完成 =====\n');
fprintf('  目录: %s\n', out_dir);
fprintf('  Step1: fig01~04 (run_scenario) + 副本至 figures/\n');
fprintf('  Step2: fig05~07\n');
fprintf('  Step3: fig08~12\n');
end

%% ══════════════════════════════════════════════════════════════════════
function [rd_xy, ac_xy, n_rd, n_ac] = extract_all_detections_xy(allDetections)
nF = length(allDetections);
n_rd = zeros(1, nF);
n_ac = zeros(1, nF);
rd_list = zeros(0, 2);
ac_list = zeros(0, 2);
for t = 1:nF
    [zr, za] = frame_detections_xy(allDetections{t});
    n_rd(t) = size(zr, 1);
    n_ac(t) = size(za, 1);
    if ~isempty(zr), rd_list = [rd_list; zr]; end %#ok<AGROW>
    if ~isempty(za), ac_list = [ac_list; za]; end %#ok<AGROW>
end
rd_xy = rd_list;
ac_xy = ac_list;
end

function [rd_xy, ac_xy] = frame_detections_xy(frame)
rd_xy = zeros(0, 2);
ac_xy = zeros(0, 2);
if isempty(frame), return; end
if iscell(frame), det_list = frame; else, det_list = num2cell(frame); end
for d = 1:length(det_list)
    det = det_list{d};
    if isobject(det)
        meas = det.Measurement(:);
        sid  = det.SensorIndex;
    elseif isstruct(det)
        meas = det.Measurement(:);
        if isfield(det, 'SensorIndex'), sid = det.SensorIndex; else, sid = 0; end
    else
        continue
    end
    if numel(meas) < 2, continue; end
    xy = meas(1:2)';
    if sid > 100
        rd_xy = [rd_xy; xy]; %#ok<AGROW>
    else
        ac_xy = [ac_xy; xy]; %#ok<AGROW>
    end
end
end

function plot_buildings_2d(buildings, Map)
hold on;
maxH = 1;
for i = 1:length(buildings)
    b = buildings(i);
    if isfield(b, 'pos')
        bx = b.pos(1); by = b.pos(2);
        bw = b.dim(1); bd = b.dim(2);
        bh = b.dim(3);
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

function rmse_ts = compute_rmse_timeseries(tracks, truth, time, cfg)
N = length(time);
rmse_ts = nan(1, N);
N_truth = length(truth);
for k = 1:N
    errs = [];
    for i = 1:length(tracks)
        if tracks(i).total_hits < cfg.mn_M, continue; end
        th = tracks(i).time_hist;
        [~, loc] = min(abs(th - k));
        if isempty(loc) || abs(th(loc) - k) > 1, continue; end
        est_xy = tracks(i).state_hist(1:2, loc);
        best_d = inf;
        for j = 1:N_truth
            raw = truth(j).pos3D;
            if k <= size(raw, 1) && ~any(isnan(raw(k, :)))
                d = norm(est_xy - raw(k, 1:2)');
                if d < best_d, best_d = d; end
            end
        end
        if best_d < 100, errs(end+1) = best_d^2; end %#ok<AGROW>
    end
    if ~isempty(errs), rmse_ts(k) = sqrt(mean(errs)); end
end
end

function copy_step1_figures(src_dir, out_dir)
names = {'fig1_3d_overview.png', 'fig2_topdown.png', ...
         'fig3_altitude.png', 'fig4_speed.png'};
map_to = {'fig01_step1_3d_overview.png', 'fig02_step1_topdown.png', ...
          'fig03_step1_altitude.png', 'fig04_step1_speed.png'};
for i = 1:length(names)
    src = fullfile(src_dir, names{i});
    if exist(src, 'file')
        copyfile(src, fullfile(out_dir, map_to{i}), 'f');
    end
end
end
