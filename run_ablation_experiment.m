%% =========================================================================
%  run_ablation_experiment.m
%  消融实验：原始融合 vs. 雷达优先融合（Radar-Priority Fusion）
%
%  科学假设：
%    H0 - 声学探测不显著影响 PDA 精度
%    H1 - 声学探测污染 R̄，拉低融合 RMSE 至声学水平
%
%  实验设计：
%    模式 A（原始融合）：cfg.acoustic_assist_only = false
%    模式 B（雷达优先）：cfg.acoustic_assist_only = true
%
%  用法：直接运行，结果保存至 ablation_results.mat
%        图表保存至 Evaluation/figures/ablation_*.png
% =========================================================================
clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'IMM'));
addpath(fullfile(this_dir, 'JPDA'));
addpath(fullfile(this_dir, 'EKF'));
addpath(fullfile(this_dir, 'TrackManagement'));
addpath(fullfile(this_dir, 'Utilities'));
addpath(fullfile(this_dir, 'Evaluation'));

%% ── 数据加载（与主程序相同路径逻辑）─────────────────────────────────────
scenario_dir = fullfile(this_dir, 'Scenario');
s2_path = fullfile(scenario_dir, 'Step2_HeteroDetections.mat');
s1_path = fullfile(scenario_dir, 'Step1_Data.mat');

if ~exist(s2_path, 'file') || ~exist(s1_path, 'file')
    error('缺少 Step1_Data.mat 或 Step2_HeteroDetections.mat，请先运行场景生成脚本。');
end

load(s2_path, 'allDetections', 'time');
data1 = load(s1_path);
truth = data1.truth;
N     = length(time);

fprintf('数据加载完成：%d 帧，%d 条真值轨迹\n', N, length(truth));

%% ── 运行两种模式 ─────────────────────────────────────────────────────────
modes = struct( ...
    'name',  {'原始融合 (Original)',    '雷达优先融合 (Radar-Priority)'}, ...
    'flag',  {false,                    true} );

results = struct();

for m_idx = 1 : 2
    fprintf('\n===== 模式 %d：%s =====\n', m_idx, modes(m_idx).name);

    cfg = config_tracker();
    cfg.dt                   = time(2) - time(1);
    cfg.acoustic_assist_only = modes(m_idx).flag;
    cfg.save_mat             = false;   % 消融实验不覆盖主结果
    cfg.debug                = false;

    [tracks_out, metrics_out, frame_log_out, extra] = ...
        run_one_mode(allDetections, truth, time, N, cfg);

    results(m_idx).name        = modes(m_idx).name;
    results(m_idx).tracks      = tracks_out;
    results(m_idx).metrics     = metrics_out;
    results(m_idx).frame_log   = frame_log_out;
    results(m_idx).extra       = extra;
    results(m_idx).cfg         = cfg;
end

%% ── 计算扩展指标 ─────────────────────────────────────────────────────────
fprintf('\n===== 扩展指标计算 =====\n');
for m_idx = 1 : 2
    results(m_idx).ext = compute_extended_metrics( ...
        results(m_idx).tracks, truth, time, N, results(m_idx).extra);
end

%% ── 打印对比表格 ─────────────────────────────────────────────────────────
print_comparison_table(results);

%% ── 生成对比图 ───────────────────────────────────────────────────────────
fig_dir = fullfile(this_dir, 'Evaluation', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
plot_ablation_comparison(results, time, fig_dir);

%% ── 保存结果 ─────────────────────────────────────────────────────────────
save_path = fullfile(this_dir, 'ablation_results.mat');
save(save_path, 'results', 'time', 'truth', '-v7.3');
fprintf('\n消融实验结果已保存：%s\n', save_path);


%% =========================================================================
%% 子函数
%% =========================================================================

function [tracks, metrics, frame_log, extra] = run_one_mode(allDetections, truth, time, N, cfg)
% 运行跟踪器主循环（与 main_imm_jpda_ekf.m 逻辑一致）

tracks  = struct([]);
next_id = 1;
frame_log = struct('n_confirmed', zeros(1,N), ...
                   'n_tentative', zeros(1,N), ...
                   'n_deleted',   zeros(1,N));

% 记录每帧传感器分布
extra.radar_only_frames   = 0;
extra.acoustic_only_frames = 0;
extra.both_frames         = 0;
extra.no_det_frames       = 0;
extra.radar_priority_activations = 0;   % 雷达优先模式触发次数（per track-frame）

use_radar_priority = cfg.acoustic_assist_only;

for k = 1 : N
    %% 步骤 1：IMM 预测
    for i = 1 : length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end
        [tracks(i).models, c_bar_i] = imm_mix(tracks(i).models, tracks(i).mu, tracks(i).Pi);
        tracks(i).c_bar_ = c_bar_i;
        tracks(i).models = imm_predict(tracks(i).models);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% 步骤 2&3：量测 + JPDA
    [Z_k, R_list_k, sensor_types_k] = detections_to_ZR(allDetections{k}, cfg);

    % 统计传感器分布
    n_radar = sum(cellfun(@(s) strcmpi(s,'Radar'), sensor_types_k));
    n_acou  = sum(cellfun(@(s) ~strcmpi(s,'Radar'), sensor_types_k));
    if isempty(sensor_types_k)
        extra.no_det_frames = extra.no_det_frames + 1;
    elseif n_radar > 0 && n_acou > 0
        extra.both_frames = extra.both_frames + 1;
    elseif n_radar > 0
        extra.radar_only_frames = extra.radar_only_frames + 1;
    else
        extra.acoustic_only_frames = extra.acoustic_only_frames + 1;
    end

    if isempty(tracks)
        active_mask = false(1, 0);
    else
        active_mask = ~strcmp({tracks.status}, 'deleted');
    end
    active_idx = find(active_mask);
    n_active   = length(active_idx);

    if n_active > 0 && size(Z_k, 2) > 0
        active_tracks = tracks(active_idx);
        assoc = jpda_run(active_tracks, Z_k, R_list_k, cfg);
    else
        assoc.valid_mat  = false(n_active, size(Z_k,2));
        assoc.innov_data = cell(n_active, size(Z_k,2));
        assoc.beta       = zeros(n_active, size(Z_k,2));
        assoc.beta0      = ones(n_active, 1);
        assoc.clusters   = {};
    end

    %% 步骤 4：PDA 更新（含雷达优先逻辑）
    is_radar_vec = cellfun(@(s) strcmpi(s,'Radar'), sensor_types_k);

    for ti = 1 : n_active
        i = active_idx(ti);

        if use_radar_priority && ~isempty(sensor_types_k)
            radar_in_gate = assoc.valid_mat(ti, :) & is_radar_vec;
            if any(radar_in_gate)
                assoc_ti = assoc;
                acoustic_mask = ~is_radar_vec;
                assoc_ti.valid_mat(ti, acoustic_mask) = false;
                assoc_ti.beta(ti,  acoustic_mask)     = 0;
                sum_beta_ti = sum(assoc_ti.beta(ti, :));
                assoc_ti.beta0(ti) = max(1 - sum_beta_ti, 0);
                extra.radar_priority_activations = extra.radar_priority_activations + 1;
            else
                assoc_ti = assoc;
            end
        else
            assoc_ti = assoc;
        end

        tracks(i) = track_imm_pda_update(tracks(i), Z_k, R_list_k, assoc_ti, ti, cfg);
    end

    %% 步骤 5-7：轨迹管理
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, k, next_id, cfg);

    %% 步骤 8：帧统计
    if ~isempty(tracks)
        statuses = {tracks.status};
        frame_log.n_confirmed(k) = sum(strcmp(statuses, 'confirmed'));
        frame_log.n_tentative(k) = sum(strcmp(statuses, 'tentative'));
        frame_log.n_deleted(k)   = sum(strcmp(statuses, 'deleted'));
    end

    if ~isempty(tracks) && isfield(tracks, 'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end
end

metrics = eval_metrics(tracks, truth, time, cfg);
end


function ext = compute_extended_metrics(tracks, truth, time, N, extra)
% 计算覆盖率、轨迹连续性、碎片化率等扩展指标

ext.radar_only_frames    = extra.radar_only_frames;
ext.acoustic_only_frames = extra.acoustic_only_frames;
ext.both_frames          = extra.both_frames;
ext.no_det_frames        = extra.no_det_frames;
ext.radar_priority_activations = extra.radar_priority_activations;

% 统计已确认轨迹
conf_mask = false(1, length(tracks));
for i = 1 : length(tracks)
    if strcmp(tracks(i).status, 'confirmed') || tracks(i).total_hits >= 3
        conf_mask(i) = true;
    end
end
conf_idx = find(conf_mask);

% 轨迹覆盖率：已确认轨迹在所有帧中"有输出"的帧数 / 真值总帧数
total_truth_frames = N * length(truth);
covered_frames = 0;
for i = conf_idx
    if ~isempty(tracks(i).time_hist)
        covered_frames = covered_frames + length(tracks(i).time_hist);
    end
end
ext.coverage_rate = min(covered_frames / max(total_truth_frames, 1), 1.0);

% 平均轨迹存活时间（帧）
survival_times = [];
for i = conf_idx
    if ~isempty(tracks(i).time_hist)
        survival_times(end+1) = length(tracks(i).time_hist); %#ok<AGROW>
    end
end
if ~isempty(survival_times)
    ext.mean_survival = mean(survival_times);
    ext.min_survival  = min(survival_times);
else
    ext.mean_survival = 0;
    ext.min_survival  = 0;
end

% 轨迹碎片化：确认后被删除且存活时间短（< 15 帧）的比例
short_tracks = 0;
for i = conf_idx
    if strcmp(tracks(i).status, 'deleted') && tracks(i).age < 15
        short_tracks = short_tracks + 1;
    end
end
ext.frag_rate = short_tracks / max(length(conf_idx), 1);

% 虚假轨迹数（确认但很快删除）
false_count = 0;
for i = conf_idx
    if strcmp(tracks(i).status, 'deleted') && tracks(i).age < 15
        false_count = false_count + 1;
    end
end
ext.false_track_count = false_count;
ext.confirmed_count   = length(conf_idx);
end


function print_comparison_table(results)
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║           消融实验对比结果                                       ║\n');
fprintf('╠══════════════════════════════╦═════════════════╦═════════════════╣\n');
fprintf('║ 指标 (Metric)                ║ 原始融合        ║ 雷达优先融合    ║\n');
fprintf('╠══════════════════════════════╬═════════════════╬═════════════════╣\n');

m_orig = results(1).metrics;
m_prio = results(2).metrics;
e_orig = results(1).ext;
e_prio = results(2).ext;

% RMSE
rmse_A = m_orig.pos_rmse; if isnan(rmse_A), rmse_A = -1; end
rmse_B = m_prio.pos_rmse; if isnan(rmse_B), rmse_B = -1; end
fprintf('║ 位置 RMSE (m)                ║ %13.2f   ║ %13.2f   ║\n', rmse_A, rmse_B);

% RMSE 变化
if rmse_A > 0 && rmse_B > 0
    delta_pct = (rmse_B - rmse_A) / rmse_A * 100;
    fprintf('║ RMSE 变化 (%%)               ║       基准       ║ %+12.1f%%  ║\n', delta_pct);
end

fprintf('╠══════════════════════════════╬═════════════════╬═════════════════╣\n');
fprintf('║ 已确认轨迹数                  ║ %13d   ║ %13d   ║\n', ...
    m_orig.confirmed_tracks, m_prio.confirmed_tracks);
fprintf('║ 虚假轨迹数                    ║ %13d   ║ %13d   ║\n', ...
    m_orig.false_tracks, m_prio.false_tracks);
fprintf('║ 平均确认延迟 (帧)             ║ %13.1f   ║ %13.1f   ║\n', ...
    nanmean_safe(m_orig.confirm_delay), nanmean_safe(m_prio.confirm_delay));
fprintf('╠══════════════════════════════╬═════════════════╬═════════════════╣\n');
fprintf('║ 轨迹覆盖率 (%%)               ║ %13.1f   ║ %13.1f   ║\n', ...
    e_orig.coverage_rate*100, e_prio.coverage_rate*100);
fprintf('║ 平均轨迹存活时间 (帧)         ║ %13.1f   ║ %13.1f   ║\n', ...
    e_orig.mean_survival, e_prio.mean_survival);
fprintf('║ 轨迹碎片化率 (%%)             ║ %13.1f   ║ %13.1f   ║\n', ...
    e_orig.frag_rate*100, e_prio.frag_rate*100);
fprintf('╠══════════════════════════════╬═════════════════╬═════════════════╣\n');
fprintf('║ 雷达优先触发次数 (track-帧)  ║      N/A         ║ %13d   ║\n', ...
    e_prio.radar_priority_activations);
fprintf('╚══════════════════════════════╩═════════════════╩═════════════════╝\n\n');

% 结论
if rmse_A > 0 && rmse_B > 0
    if rmse_B < rmse_A * 0.5
        fprintf('[结论] H1 获得支持：声学探测正在主导 PDA 更新，雷达优先融合 RMSE 显著下降 (%.1f%%)。\n', ...
            (rmse_A-rmse_B)/rmse_A*100);
        fprintf('       建议保留 Radar-Priority Fusion 模式作为默认策略。\n');
    elseif rmse_B < rmse_A * 0.9
        fprintf('[结论] H1 部分支持：声学探测对精度有一定影响，雷达优先融合有所改善 (%.1f%%)。\n', ...
            (rmse_A-rmse_B)/rmse_A*100);
        fprintf('       需结合覆盖率权衡是否启用。\n');
    else
        fprintf('[结论] H0 获得支持：RMSE 无显著差异，融合 RMSE 退化另有原因。\n');
        fprintf('       真正瓶颈分析：见下方诊断输出。\n');
        diagnose_rmse_bottleneck(results);
    end
end
end


function diagnose_rmse_bottleneck(results)
fprintf('\n===== RMSE 瓶颈诊断 =====\n');
e = results(1).ext;
fprintf('  传感器覆盖统计：\n');
fprintf('    仅雷达帧：    %d\n', e.radar_only_frames);
fprintf('    仅声学帧：    %d\n', e.acoustic_only_frames);
fprintf('    雷达+声学帧： %d\n', e.both_frames);
fprintf('    无探测帧：    %d\n', e.no_det_frames);
fprintf('\n  可能原因分析：\n');
if e.acoustic_only_frames > e.radar_only_frames
    fprintf('  → 声学专属帧数较多，轨迹在这些帧被声学拉偏后雷达帧累积误差\n');
end
if e.both_frames > 0
    fprintf('  → 雷达+声学共存帧 %d 个，声学在这些帧污染 R̄\n', e.both_frames);
end
fprintf('  → 建议检查：eval_metrics 中距离上限 100m 是否过宽，导致错配轨迹计入 RMSE\n');
end


function plot_ablation_comparison(results, time, fig_dir)
N = length(time);

figure('Position', [100 100 1200 800], 'Visible', 'off');

% 子图 1：逐帧已确认轨迹数
subplot(2, 3, 1);
plot(time, results(1).frame_log.n_confirmed, 'b-', 'LineWidth', 1.5); hold on;
plot(time, results(2).frame_log.n_confirmed, 'r--', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('已确认轨迹数');
title('轨迹确认数对比');
legend('原始融合', '雷达优先融合', 'Location', 'best');
grid on;

% 子图 2：轨迹存活时间分布（直方图）
subplot(2, 3, 2);
sv_A = get_survival_times(results(1).tracks);
sv_B = get_survival_times(results(2).tracks);
if ~isempty(sv_A) || ~isempty(sv_B)
    edges = 0:5:max([sv_A(:); sv_B(:); 1])+5;
    histogram(sv_A, edges, 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
    histogram(sv_B, edges, 'FaceColor', 'r', 'FaceAlpha', 0.5);
end
xlabel('存活时间 (帧)'); ylabel('轨迹数量');
title('轨迹存活时间分布');
legend('原始融合', '雷达优先融合', 'Location', 'best');
grid on;

% 子图 3：RMSE 对比柱状图
subplot(2, 3, 3);
rmse_vals = [results(1).metrics.pos_rmse, results(2).metrics.pos_rmse];
bar_h = bar(rmse_vals, 'FaceColor', 'flat');
bar_h.CData = [0 0.4 0.8; 0.8 0.2 0.2];
set(gca, 'XTickLabel', {'原始融合', '雷达优先'});
ylabel('位置 RMSE (m)');
title('RMSE 对比');
grid on;
for bi = 1:2
    text(bi, rmse_vals(bi) + 0.1, sprintf('%.2f m', rmse_vals(bi)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

% 子图 4：覆盖率对比
subplot(2, 3, 4);
cov_vals = [results(1).ext.coverage_rate, results(2).ext.coverage_rate] * 100;
bar_h2 = bar(cov_vals, 'FaceColor', 'flat');
bar_h2.CData = [0 0.4 0.8; 0.8 0.2 0.2];
set(gca, 'XTickLabel', {'原始融合', '雷达优先'});
ylabel('覆盖率 (%)'); ylim([0 105]);
title('轨迹覆盖率对比');
grid on;

% 子图 5：虚假轨迹 & 碎片化
subplot(2, 3, 5);
fa_vals   = [results(1).metrics.false_tracks, results(2).metrics.false_tracks];
frag_vals = [results(1).ext.frag_rate, results(2).ext.frag_rate] * 100;
yyaxis left;
bar(fa_vals, 'FaceAlpha', 0.6);
ylabel('虚假轨迹数');
yyaxis right;
plot(1:2, frag_vals, 'ko-', 'MarkerFaceColor', 'k');
ylabel('碎片化率 (%)');
set(gca, 'XTickLabel', {'原始融合', '雷达优先'});
title('虚假轨迹 & 碎片化');
grid on;

% 子图 6：传感器帧类型分布（饼图，原始模式）
subplot(2, 3, 6);
e = results(1).ext;
pie_vals = [e.radar_only_frames, e.acoustic_only_frames, e.both_frames, e.no_det_frames];
pie_vals = max(pie_vals, 0);
if sum(pie_vals) > 0
    pie(pie_vals + 1e-9, {'仅雷达', '仅声学', '雷达+声学', '无探测'});
end
title('传感器帧类型分布');

sgtitle('消融实验对比：原始融合 vs. 雷达优先融合', 'FontSize', 14, 'FontWeight', 'bold');

% 保存
fname = fullfile(fig_dir, 'ablation_comparison.png');
print(gcf, fname, '-dpng', '-r300');
fprintf('对比图已保存：%s\n', fname);
close(gcf);
end


function sv = get_survival_times(tracks)
sv = [];
for i = 1 : length(tracks)
    if (strcmp(tracks(i).status, 'confirmed') || tracks(i).total_hits >= 3) && ...
       ~isempty(tracks(i).time_hist)
        sv(end+1) = length(tracks(i).time_hist); %#ok<AGROW>
    end
end
end


function v = nanmean_safe(x)
if isnan(x) || isempty(x)
    v = NaN;
else
    v = x;
end
end
