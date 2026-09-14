%% run_acoustic_likelihood_ab.m
% PF-JPDA-IMM 声学似然 A/B 验证（post vs L_softmax）
% 不修改任何算法源码；A 组通过 validation/pf_jpda_update_ab.m 复现旧行为。
%
% 用法（MATLAB 命令行）：
%   cd Tracker
%   run('validation/run_acoustic_likelihood_ab.m')
%
% 产出：
%   validation/ab_results.mat
%   report/PF_acoustic_likelihood_AB_验证报告.md
%   report/figures/ab_*.png

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
tracker_dir = fileparts(this_dir);
addpath(tracker_dir);
addpath(fullfile(tracker_dir, 'IMM'));
addpath(fullfile(tracker_dir, 'JPDA'));
addpath(fullfile(tracker_dir, 'EKF'));
addpath(fullfile(tracker_dir, 'PF'));
addpath(fullfile(tracker_dir, 'TrackManagement'));
addpath(fullfile(tracker_dir, 'Utilities'));
addpath(fullfile(tracker_dir, 'Evaluation'));
addpath(fullfile(tracker_dir, 'Scenario'));
addpath(fullfile(tracker_dir, 'matlab'));
addpath(this_dir);

fig_dir = fullfile(tracker_dir, 'report', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fprintf('===== PF 声学似然 A/B 验证 =====\n');
fprintf('A = ac_pack.post（旧版）\n');
fprintf('B = ac_pack.L_softmax（当前生产代码 pf_jpda_update.m）\n\n');

%% ── 数据加载 ─────────────────────────────────────────────────────────────
[s1, s2, data1, s2data] = load_scenario_data(tracker_dir);
truth         = data1.truth;
allDetections = s2data.allDetections;
time          = s2data.time;
N             = length(time);

cfg = config_tracker();
cfg.use_pf   = true;
cfg.save_mat = false;
cfg.dt       = time(2) - time(1);
ac_cfg       = acoustic_config();

fprintf('  帧数: %d  目标: %d  N_particle: %d  resample_thresh: %.2f\n\n', ...
    N, length(truth), cfg.pf_N, cfg.pf_resample_thresh);

%% ── 代码确认（静态） ─────────────────────────────────────────────────────
code_check = verify_code_modification(tracker_dir);
print_code_check(code_check);

%% ── A/B 运行 ─────────────────────────────────────────────────────────────
variants(1) = struct('label', 'A_post', ...
    'title', 'A (post)', ...
    'use_prod_pf', false, ...
    'ac_field_mode', 'post');
variants(2) = struct('label', 'B_L_softmax', ...
    'title', 'B (L_softmax)', ...
    'use_prod_pf', true, ...
    'ac_field_mode', 'L_softmax');

% 用 cell 存两次运行结果，避免 tracks 字段不一致导致 struct 数组赋值失败
results = cell(1, 2);
for vi = 1:2
    fprintf('\n===== 运行 %s =====\n', variants(vi).title);
    rng(42, 'twister');
    results{vi} = run_one_variant( ...
        allDetections, truth, time, N, cfg, ac_cfg, variants(vi));
end

%% ── 汇总指标 ─────────────────────────────────────────────────────────────
summary = build_summary(results, cfg, N);
print_summary_table(summary);

%% ── 绘图 ─────────────────────────────────────────────────────────────────
plot_ab_figures(results, summary, time, cfg, fig_dir);

%% ── 生成 Markdown 报告 ───────────────────────────────────────────────────
report_path = fullfile(tracker_dir, 'report', 'PF_acoustic_likelihood_AB_验证报告.md');
write_markdown_report(report_path, code_check, results, summary, cfg, N, fig_dir);
fprintf('\n验证报告已写入: %s\n', report_path);

save(fullfile(this_dir, 'ab_results.mat'), ...
    'results', 'summary', 'code_check', 'cfg', 'time', 'truth', '-v7.3');
fprintf('原始结果已保存: validation/ab_results.mat\n');

%% ── 最终建议（控制台） ───────────────────────────────────────────────────
rec = decide_recommendation(summary);
fprintf('\n【最终建议】%s\n', rec.text);


%% =========================================================================
function [s1, s2, data1, s2data] = load_scenario_data(tracker_dir)
scenario_dir = fullfile(tracker_dir, 'Scenario');
data_dir_v3  = fullfile(tracker_dir, '..', 'v3');
data_dir_v2  = fullfile(tracker_dir, '..', 'v2.0');

s2_cands = { fullfile(scenario_dir, 'Step2_HeteroDetections.mat'), ...
             fullfile(data_dir_v3,   'Step2_HeteroDetections_v3.mat'), ...
             fullfile(data_dir_v2,   'Step2_HeteroDetections.mat') };
s1_cands = { fullfile(scenario_dir, 'Step1_Data.mat'), ...
             fullfile(data_dir_v3,   'Step1_Data_v3.mat'), ...
             fullfile(data_dir_v2,   'Step1_Data.mat') };

s2 = ''; s1 = '';
for fi = 1:length(s2_cands)
    if exist(s2_cands{fi}, 'file')
        s2 = s2_cands{fi}; break;
    end
end
for fi = 1:length(s1_cands)
    if exist(s1_cands{fi}, 'file')
        s1 = s1_cands{fi}; break;
    end
end
if isempty(s2) || isempty(s1)
    error(['未找到 Step1/Step2 数据。请先运行 Scenario/run_scenario.m ', ...
           '与 step2_detection_simulation.m']);
end
s2data = load(s2, 'allDetections', 'time');
data1  = load(s1);
end


function code_check = verify_code_modification(tracker_dir)
code_check.pf_main_uses_L_softmax = false;
code_check.pf_debug_uses_L_softmax = false;
code_check.other_pf_paths_use_post = {};
code_check.call_chain_ok = false;

pf_main = fileread(fullfile(tracker_dir, 'PF', 'pf_jpda_update.m'));
pf_dbg  = fileread(fullfile(tracker_dir, 'PF', 'pf_jpda_update_debug.m'));

code_check.pf_main_uses_L_softmax  = contains(pf_main, 'L_ac_field = ac_pack.L_softmax');
code_check.pf_debug_uses_L_softmax = contains(pf_dbg,  'L_ac_field = ac_pack.L_softmax');
code_check.pf_main_no_post_field   = ~contains(pf_main, 'L_ac_field = ac_pack.post');
code_check.pf_debug_no_post_field  = ~contains(pf_dbg,  'L_ac_field = ac_pack.post');

% 扫描 PF 更新路径
scan_files = { ...
    fullfile(tracker_dir, 'PF', 'pf_jpda_update.m'), ...
    fullfile(tracker_dir, 'PF', 'pf_jpda_update_debug.m'), ...
    fullfile(tracker_dir, 'TrackManagement', 'track_imm_pda_update.m'), ...
    fullfile(tracker_dir, 'TrackManagement', 'track_imm_pf_update.m') };
for fi = 1:length(scan_files)
    txt = fileread(scan_files{fi});
    if contains(txt, 'ac_pack.post') && contains(txt, 'L_ac_field')
        code_check.other_pf_paths_use_post{end+1} = scan_files{fi}; %#ok<AGROW>
    end
end

ac_txt = fileread(fullfile(tracker_dir, 'Scenario', 'acoustic_field_to_tracker.m'));
code_check.acoustic_to_tracker_exists = exist(fullfile(tracker_dir, 'Scenario', 'acoustic_field_to_tracker.m'), 'file') == 2;
code_check.ac_pack_has_post_for_recursion = contains(ac_txt, 'prev_post_out = ac_pack.post');

code_check.call_chain_ok = code_check.acoustic_to_tracker_exists && ...
    code_check.pf_main_uses_L_softmax && code_check.pf_debug_uses_L_softmax && ...
    code_check.pf_main_no_post_field && isempty(code_check.other_pf_paths_use_post);
end


function print_code_check(c)
fprintf('--- 第一部分：修改确认 ---\n');
fprintf('  pf_jpda_update.m 使用 L_softmax: %s\n', bool_str(c.pf_main_uses_L_softmax));
fprintf('  pf_jpda_update_debug.m 使用 L_softmax: %s\n', bool_str(c.pf_debug_uses_L_softmax));
fprintf('  其他 PF 路径仍用 ac_pack.post 作似然: %s\n', ...
    iif(isempty(c.other_pf_paths_use_post), '否', strjoin(c.other_pf_paths_use_post, '; ')));
fprintf('  acoustic_field_to_tracker → ac_pack → pf_jpda_update 链: %s\n', bool_str(c.call_chain_ok));
fprintf('  ac_pack.post 仍用于帧间递归 (prev_post): %s\n\n', bool_str(c.ac_pack_has_post_for_recursion));
end


function out = run_one_variant(allDetections, truth, time, N, cfg, ac_cfg, variant)
tracks   = struct([]);
next_id  = 1;
prev_post_ac = [];

log = init_run_log(N, cfg);

snap_k = pick_snapshot_frames(N, truth, allDetections);
log.snap_frames = snap_k;

frame_log = struct('n_confirmed', zeros(1,N), ...
                   'n_tentative', zeros(1,N), ...
                   'n_deleted',   zeros(1,N));

for k = 1:N
    if mod(k, 100) == 0
        fprintf('  帧 %4d/%d\n', k, N);
    end

    ac_pack_k = [];
    src_pos_k = [];
    for tn = 1:length(truth)
        if k <= size(truth(tn).pos3D, 1) && ~any(isnan(truth(tn).pos3D(k,:)))
            src_pos_k = [src_pos_k; truth(tn).pos3D(k, 1:3)]; %#ok<AGROW>
        end
    end
    if ~isempty(src_pos_k)
        [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
            src_pos_k, ac_cfg, prev_post_ac, 'softmax');
    end

    for i = 1:length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end
        [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
            tracks(i).mu, tracks(i).Pi);
        tracks(i).c_bar_ = c_bar_i;
        tracks(i).models = imm_pf_predict(tracks(i).models);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    [Z_k, R_list_k, sensor_types_k] = detections_to_ZR(allDetections{k}, cfg);
    D = size(Z_k, 2);

    if isempty(tracks)
        active_mask = false(1, 0);
    else
        active_mask = ~strcmp({tracks.status}, 'deleted');
    end
    active_idx = find(active_mask);
    n_active   = length(active_idx);

    if n_active > 0 && D > 0
        assoc = jpda_run(tracks(active_idx), Z_k, R_list_k, cfg);
    else
        assoc.valid_mat  = false(n_active, D);
        assoc.innov_data = cell(n_active, D);
        assoc.beta       = zeros(n_active, D);
        assoc.beta0      = ones(n_active, 1);
        assoc.clusters   = {};
    end

    log = update_jpda_log(log, k, assoc, n_active, D, sensor_types_k);

    frame_ess = [];
    frame_mu1 = [];
    frame_mu2 = [];
    frame_rs  = 0;

    for ti = 1:n_active
        i = active_idx(ti);
        [tracks(i), fs] = track_imm_pda_update_ab( ...
            tracks(i), Z_k, R_list_k, assoc, ti, cfg, ac_pack_k, ...
            variant.use_prod_pf, variant.ac_field_mode);

        frame_ess(end+1) = fs.ess; %#ok<AGROW>
        frame_mu1(end+1) = fs.mu(1); %#ok<AGROW>
        frame_mu2(end+1) = fs.mu(2); %#ok<AGROW>
        frame_rs = frame_rs + fs.n_resample;

        if ismember(k, snap_k)
            si = find(snap_k == k, 1);
            log.snap(si).frame = k;
            nt = struct('id', tracks(i).id, ...
                        'particles_xy', fs.particles_xy, ...
                        'weights', fs.weights, ...
                        'mu', fs.mu(:)', ...
                        'ess', fs.ess, ...
                        'dominant_model', fs.dominant_model);
            log.snap(si).tracks = [log.snap(si).tracks, nt]; %#ok<AGROW>
        end
    end

    if ~isempty(frame_ess)
        log.ess_k(k)       = min(frame_ess);
        log.ess_mean_k(k)  = mean(frame_ess);
        log.mu1_k(k)       = mean(frame_mu1);
        log.mu2_k(k)       = mean(frame_mu2);
        log.entropy_k(k)   = imm_entropy(mean([frame_mu1; frame_mu2], 2));
        log.resample_k(k)  = frame_rs;
        if any(log.resample_k(max(1,k-1):k) > 0) && log.last_resample_frame > 0
            log.resample_intervals(end+1) = k - log.last_resample_frame; %#ok<AGROW>
        end
        if frame_rs > 0
            log.last_resample_frame = k;
            log.total_resamples = log.total_resamples + frame_rs;
        end
        log.active_update_frames = log.active_update_frames + 1;
    end

    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, k, next_id, cfg);

    if ~isempty(tracks)
        st = {tracks.status};
        frame_log.n_confirmed(k) = sum(strcmp(st,'confirmed'));
        frame_log.n_tentative(k) = sum(strcmp(st,'tentative'));
        frame_log.n_deleted(k)   = sum(strcmp(st,'deleted'));
    end

    if ~isempty(tracks) && isfield(tracks, 'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end
end

metrics = eval_metrics(tracks, truth, time, cfg);
log = finalize_run_log(log, tracks, cfg, metrics, frame_log);

out.label       = variant.label;
out.title       = variant.title;
out.metrics     = metrics;
out.frame_log   = frame_log;
out.log         = log;
% 仅保留轨迹摘要，完整 tracks 各 run 字段集可能不同，不宜放入 struct 数组
out.track_summary = summarize_tracks_for_ab(tracks, metrics, frame_log);
end


function log = init_run_log(N, cfg)
log.ess_k              = nan(N, 1);
log.ess_mean_k         = nan(N, 1);
log.ess_all            = [];
log.mu1_k              = nan(N, 1);
log.mu2_k              = nan(N, 1);
log.entropy_k          = nan(N, 1);
log.resample_k         = zeros(N, 1);
log.gate_meas_mean_k   = nan(N, 1);
log.beta0_mean_k       = nan(N, 1);
log.beta_d_max_mean_k  = nan(N, 1);
log.gate_reject_rate_k = nan(N, 1);
log.assoc_fail_k       = zeros(N, 1);
log.total_resamples    = 0;
log.resample_intervals = [];
log.last_resample_frame = 0;
log.active_update_frames = 0;
log.model_switch_count = 0;
log.snap_frames = [];
empty_snap_track = struct('id', {}, 'particles_xy', {}, 'weights', {}, ...
    'mu', {}, 'ess', {}, 'dominant_model', {});
log.snap = repmat(struct('frame', [], 'tracks', empty_snap_track), 1, 4);
log.cfg_N = cfg.pf_N;
log.cfg_thresh = cfg.pf_resample_thresh;
end


function snap_k = pick_snapshot_frames(N, truth, allDetections)
% 起始 / 稳态 / 机动 / 声学峰值 四个关键时刻
k_start = max(1, round(0.1 * N));
k_mid   = round(0.5 * N);
k_late  = round(0.75 * N);

% 声学峰值：acoustic 探测数最多帧
ac_counts = zeros(1, numel(allDetections));
for k = 1:numel(allDetections)
    frame = allDetections{k};
    if isempty(frame), continue; end
    for d = 1:length(frame)
        if isfield(frame(d), 'sensorType') && ~strcmpi(frame(d).sensorType, 'Radar')
            ac_counts(k) = ac_counts(k) + 1;
        end
    end
end
[~, k_ac_peak] = max(ac_counts);
if k_ac_peak < 1, k_ac_peak = k_mid; end

snap_k = unique([k_start, k_mid, k_late, k_ac_peak]);
end


function log = update_jpda_log(log, k, assoc, n_active, D, sensor_types_k)
if n_active == 0 || D == 0
    log.gate_reject_rate_k(k) = NaN;
    return;
end

gate_counts = sum(assoc.valid_mat, 2);
log.gate_meas_mean_k(k)  = mean(gate_counts);
log.beta0_mean_k(k)      = mean(assoc.beta0);
beta_max = zeros(n_active, 1);
for ti = 1:n_active
    row = assoc.beta(ti, :);
    if any(row > 0)
        beta_max(ti) = max(row);
    end
end
log.beta_d_max_mean_k(k) = mean(beta_max);

total_pairs = n_active * D;
in_gate     = sum(assoc.valid_mat(:));
log.gate_reject_rate_k(k) = 1 - in_gate / max(total_pairs, 1);

% 关联失败：有量测但所有轨迹 beta 为 0 且 beta0≈1
if D > 0 && all(sum(assoc.beta, 2) < 1e-12)
    log.assoc_fail_k(k) = 1;
end
end


function H = imm_entropy(mu)
mu = max(mu(:), 1e-300);
mu = mu / sum(mu);
H = -sum(mu .* log(mu));
end


function log = finalize_run_log(log, tracks, cfg, metrics, frame_log)
ess_all = [];
mu_switch = 0;
for i = 1:length(tracks)
    if ~isempty(tracks(i).ess_hist)
        ess_all = [ess_all, tracks(i).ess_hist(:)']; %#ok<AGROW>
    end
    if isfield(tracks(i), 'mu_hist') && ~isempty(tracks(i).mu_hist)
        mh = tracks(i).mu_hist;
        for t = 2:size(mh, 2)
            [~, p1] = max(mh(:, t-1));
            [~, p2] = max(mh(:, t));
            if p1 ~= p2, mu_switch = mu_switch + 1; end
        end
    end
end
log.ess_all = ess_all;

if ~isempty(ess_all)
    log.ess_mean  = mean(ess_all);
    log.ess_min   = min(ess_all);
    log.ess_std   = std(ess_all);
    log.ess_below = 100 * mean(ess_all < cfg.pf_resample_thresh * cfg.pf_N);
else
    log.ess_mean = NaN; log.ess_min = NaN; log.ess_std = NaN; log.ess_below = NaN;
end

Nf = log.active_update_frames;
log.resample_ratio = log.total_resamples / max(Nf * cfg.n_models, 1);
if ~isempty(log.resample_intervals)
    log.resample_interval_mean = mean(log.resample_intervals);
    log.resample_interval_min  = min(log.resample_intervals);
else
    log.resample_interval_mean = NaN;
    log.resample_interval_min  = NaN;
end

log.mu1_mean = nanmean(log.mu1_k);
log.mu2_mean = nanmean(log.mu2_k);
log.entropy_mean = nanmean(log.entropy_k);
log.model_switch_count = mu_switch;

log.metrics = metrics;
log.frame_log = frame_log;
log.track_continuity = compute_track_continuity(tracks, frame_log);
log.lost_track_count = count_lost_tracks(tracks);
log.confirmed_count  = metrics.confirmed_tracks;
log.deleted_count    = sum(strcmp({tracks.status}, 'deleted'));
end


function ts = summarize_tracks_for_ab(tracks, metrics, frame_log)
ts.confirmed_count  = metrics.confirmed_tracks;
ts.false_count      = metrics.false_tracks;
ts.deleted_count    = sum(strcmp({tracks.status}, 'deleted'));
ts.lost_track_count = count_lost_tracks(tracks);
ts.track_continuity = compute_track_continuity(tracks, frame_log);
ts.total_tracks     = length(tracks);
end


function cont = compute_track_continuity(tracks, ~)
% 已确认轨迹 time_hist 无中断比例
ratios = [];
for i = 1:length(tracks)
    if strcmp(tracks(i).status, 'deleted') && tracks(i).total_hits < 3
        continue;
    end
    th = tracks(i).time_hist;
    if length(th) < 2, continue; end
    gaps = diff(th);
    ratios(end+1) = mean(gaps <= 1.5); %#ok<AGROW>
end
if isempty(ratios)
    cont = NaN;
else
    cont = mean(ratios);
end
end


function n = count_lost_tracks(tracks)
n = 0;
for i = 1:length(tracks)
    if strcmp(tracks(i).status, 'deleted') && tracks(i).total_hits >= 3
        n = n + 1;
    end
end
end


function summary = build_summary(results, cfg, N)
for vi = 1:2
    R = results{vi};
    L = R.log;
    M = R.metrics;
    s(vi).label = R.label; %#ok<AGROW>
    s(vi).title = R.title;
    s(vi).ess_mean = L.ess_mean;
    s(vi).ess_min = L.ess_min;
    s(vi).ess_std = L.ess_std;
    s(vi).ess_below_pct = L.ess_below;
    s(vi).total_resamples = L.total_resamples;
    s(vi).resample_ratio = L.resample_ratio;
    s(vi).resample_interval_mean = L.resample_interval_mean;
    s(vi).resample_interval_min = L.resample_interval_min;
    s(vi).mu1_mean = L.mu1_mean;
    s(vi).mu2_mean = L.mu2_mean;
    s(vi).entropy_mean = L.entropy_mean;
    s(vi).model_switch_count = L.model_switch_count;
    s(vi).pos_rmse = M.pos_rmse;
    s(vi).vel_rmse = M.vel_rmse;
    s(vi).max_pos_err = max_rmse_frame(M.rmse_pos_k);
    s(vi).mean_abs_err = nanmean(M.rmse_pos_k);
    s(vi).gate_meas_mean = nanmean(L.gate_meas_mean_k);
    s(vi).beta0_mean = nanmean(L.beta0_mean_k);
    s(vi).beta_d_max_mean = nanmean(L.beta_d_max_mean_k);
    s(vi).assoc_fail_count = sum(L.assoc_fail_k);
    s(vi).gate_reject_rate = nanmean(L.gate_reject_rate_k);
    s(vi).track_continuity = L.track_continuity;
    s(vi).lost_track_count = L.lost_track_count;
    s(vi).confirmed_count = L.confirmed_count;
    s(vi).deleted_count = L.deleted_count;
    s(vi).ess_k = L.ess_k;
    s(vi).mu1_k = L.mu1_k;
    s(vi).mu2_k = L.mu2_k;
    s(vi).entropy_k = L.entropy_k;
end
summary = s;
summary.cfg = cfg;
summary.N = N;
end


function v = max_rmse_frame(rmse_k)
v = max(rmse_k(~isnan(rmse_k)));
if isempty(v), v = NaN; end
end


function print_summary_table(summary)
fprintf('\n--- ESS 对比 ---\n');
fprintf('%-22s %12s %12s\n', '指标', 'A (post)', 'B (L_softmax)');
fprintf('%-22s %12.1f %12.1f\n', '平均 ESS', summary(1).ess_mean, summary(2).ess_mean);
fprintf('%-22s %12.1f %12.1f\n', '最小 ESS', summary(1).ess_min, summary(2).ess_min);
fprintf('%-22s %12.1f %12.1f\n', 'ESS 标准差', summary(1).ess_std, summary(2).ess_std);
fprintf('%-22s %11.2f%% %11.2f%%\n', 'ESS<0.5N 比例', summary(1).ess_below_pct, summary(2).ess_below_pct);

fprintf('\n--- 重采样对比 ---\n');
fprintf('%-22s %12d %12d\n', '总重采样次数', summary(1).total_resamples, summary(2).total_resamples);
fprintf('%-22s %12.4f %12.4f\n', '触发比例', summary(1).resample_ratio, summary(2).resample_ratio);
fprintf('%-22s %12.1f %12.1f\n', '平均间隔(帧)', summary(1).resample_interval_mean, summary(2).resample_interval_mean);
fprintf('%-22s %12d %12d\n', '最短间隔(帧)', summary(1).resample_interval_min, summary(2).resample_interval_min);

fprintf('\n--- 跟踪精度 ---\n');
fprintf('%-22s %12.2f %12.2f\n', '位置 RMSE (m)', summary(1).pos_rmse, summary(2).pos_rmse);
fprintf('%-22s %12.2f %12.2f\n', '逐帧 RMSE 均值', summary(1).mean_abs_err, summary(2).mean_abs_err);
fprintf('%-22s %12.2f %12.2f\n', '最大逐帧 RMSE', summary(1).max_pos_err, summary(2).max_pos_err);
end


function plot_ab_figures(results, summary, time, cfg, fig_dir)
N = length(time);
colors = [0.85 0.33 0.1; 0 0.45 0.74];

fig1 = figure('Visible', 'off', 'Position', [50 50 1200 500]);
for vi = 1:2
    subplot(1, 2, vi);
    plot(time, results{vi}.log.ess_k, 'Color', colors(vi,:), 'LineWidth', 1.2);
    yline(cfg.pf_resample_thresh * cfg.pf_N, 'k--', '0.5N');
    xlabel('时间 (s)'); ylabel('ESS');
    title(sprintf('%s  ESS', results{vi}.title));
    grid on;
end
sgtitle('ESS 时间曲线 A/B');
saveas(fig1, fullfile(fig_dir, 'ab_ess_curve.png'));
close(fig1);

fig2 = figure('Visible', 'off', 'Position', [50 50 1200 500]);
for vi = 1:2
    subplot(1, 2, vi);
    plot(time, results{vi}.log.mu1_k, '-', 'Color', colors(1,:), 'LineWidth', 1); hold on;
    plot(time, results{vi}.log.mu2_k, '-', 'Color', colors(2,:), 'LineWidth', 1);
    xlabel('时间 (s)'); ylabel('\mu');
    legend('\mu_1 (CV)', '\mu_2 (CA)', 'Location', 'best');
    title(sprintf('%s  IMM 模型概率', results{vi}.title));
    grid on;
end
sgtitle('IMM 模型概率 A/B');
saveas(fig2, fullfile(fig_dir, 'ab_imm_mu_curve.png'));
close(fig2);

fig3 = figure('Visible', 'off', 'Position', [50 50 1400 900]);
for si = 1:4
    for vi = 1:2
        subplot(4, 2, (si-1)*2 + vi);
        snap = results{vi}.log.snap(min(si, numel(results{vi}.log.snap)));
        if isempty(snap.tracks) || ~isfield(snap.tracks, 'particles_xy') || ...
                isempty(snap.tracks(1).particles_xy)
            title(sprintf('%s 帧%d (无轨迹)', results{vi}.title, snap.frame));
            continue;
        end
        tr = snap.tracks(1);
        scatter(tr.particles_xy(1,:), tr.particles_xy(2,:), 8, tr.weights, 'filled');
        colorbar; axis equal; grid on;
        title(sprintf('%s 帧%d ESS=%.0f', results{vi}.title, snap.frame, tr.ess));
    end
end
sgtitle('关键时刻粒子云（颜色=权重）');
saveas(fig3, fullfile(fig_dir, 'ab_particle_snapshots.png'));
close(fig3);

fig4 = figure('Visible', 'off', 'Position', [50 50 900 400]);
names = {summary(1).title, summary(2).title};
bar_data = [summary(1).ess_mean, summary(2).ess_mean; ...
            summary(1).total_resamples, summary(2).total_resamples; ...
            summary(1).pos_rmse, summary(2).pos_rmse];
bar(bar_data');
set(gca, 'XTickLabel', names);
legend('平均 ESS', '重采样次数', '位置 RMSE', 'Location', 'best');
title('A/B 核心指标');
grid on;
saveas(fig4, fullfile(fig_dir, 'ab_summary_bar.png'));
close(fig4);
end


function write_markdown_report(path, code_check, results, summary, cfg, N, fig_dir)
fid = fopen(path, 'w', 'n', 'UTF-8');

fprintf(fid, '# PF-JPDA-IMM 声学似然修复 A/B 验证报告\n\n');
fprintf(fid, '**生成时间：** %s  \n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '**仿真帧数：** %d  **粒子数 N：** %d  **重采样阈值：** %.2fN\n\n', ...
    N, cfg.pf_N, cfg.pf_resample_thresh);

fprintf(fid, '## 1 修改内容确认\n\n');
fprintf(fid, '| 检查项 | 结果 |\n|--------|------|\n');
fprintf(fid, '| `pf_jpda_update.m` 使用 `L_softmax` | %s |\n', bool_str(code_check.pf_main_uses_L_softmax));
fprintf(fid, '| `pf_jpda_update_debug.m` 与主版本一致 | %s |\n', ...
    bool_str(code_check.pf_debug_uses_L_softmax && code_check.pf_debug_no_post_field));
fprintf(fid, '| 其他 PF 路径仍用 `post` 作似然 | %s |\n', ...
    iif(isempty(code_check.other_pf_paths_use_post), '否', '是'));
fprintf(fid, '| 调用链 acoustic_field_to_tracker → ac_pack → pf_jpda_update | %s |\n', ...
    bool_str(code_check.call_chain_ok));
fprintf(fid, '| `post` 仍用于帧间递归 prev_post | %s |\n\n', ...
    bool_str(code_check.ac_pack_has_post_for_recursion));

fprintf(fid, '调用链：\n\n');
fprintf(fid, '```\nmain_imm_jpda_pf.m\n  └ acoustic_field_to_tracker → build_acoustic_prob_field → ac_pack\n  └ track_imm_pda_update → pf_jpda_update (L_softmax 插值)\n```\n\n');

fprintf(fid, '## 2 ESS 分析\n\n');
fprintf(fid, '| 指标 | A (post) | B (L_softmax) | Δ(B−A) |\n');
fprintf(fid, '|------|----------|---------------|--------|\n');
fprintf(fid, '| 平均 ESS | %.1f | %.1f | %+.1f |\n', ...
    summary(1).ess_mean, summary(2).ess_mean, summary(2).ess_mean - summary(1).ess_mean);
fprintf(fid, '| 最小 ESS | %.1f | %.1f | %+.1f |\n', ...
    summary(1).ess_min, summary(2).ess_min, summary(2).ess_min - summary(1).ess_min);
fprintf(fid, '| ESS 标准差 | %.1f | %.1f | %+.1f |\n', ...
    summary(1).ess_std, summary(2).ess_std, summary(2).ess_std - summary(1).ess_std);
fprintf(fid, '| ESS<0.5N 比例 (%%) | %.2f | %.2f | %+.2f |\n\n', ...
    summary(1).ess_below_pct, summary(2).ess_below_pct, ...
    summary(2).ess_below_pct - summary(1).ess_below_pct);
fprintf(fid, '![ESS曲线](../figures/ab_ess_curve.png)\n\n');

ess_conc = ess_conclusion(summary);
fprintf(fid, '**结论：** %s\n\n', ess_conc);

fprintf(fid, '## 3 IMM 分析\n\n');
fprintf(fid, '| 指标 | A (post) | B (L_softmax) |\n|------|----------|---------------|\n');
fprintf(fid, '| 平均 μ₁ | %.3f | %.3f |\n', summary(1).mu1_mean, summary(2).mu1_mean);
fprintf(fid, '| 平均 μ₂ | %.3f | %.3f |\n', summary(1).mu2_mean, summary(2).mu2_mean);
fprintf(fid, '| 模型切换次数 | %d | %d |\n', summary(1).model_switch_count, summary(2).model_switch_count);
fprintf(fid, '| 平均模型熵 H | %.3f | %.3f |\n\n', summary(1).entropy_mean, summary(2).entropy_mean);
fprintf(fid, '![IMM模型概率](../figures/ab_imm_mu_curve.png)\n\n');
fprintf(fid, '**结论：** %s\n\n', imm_conclusion(summary));

fprintf(fid, '## 4 跟踪精度分析\n\n');
fprintf(fid, '| 指标 | A (post) | B (L_softmax) |\n|------|----------|---------------|\n');
fprintf(fid, '| 位置 RMSE (m) | %.2f | %.2f |\n', summary(1).pos_rmse, summary(2).pos_rmse);
fprintf(fid, '| 逐帧 RMSE 均值 (m) | %.2f | %.2f |\n', summary(1).mean_abs_err, summary(2).mean_abs_err);
fprintf(fid, '| 最大逐帧 RMSE (m) | %.2f | %.2f |\n', summary(1).max_pos_err, summary(2).max_pos_err);
fprintf(fid, '| 轨迹连续率 | %.3f | %.3f |\n', summary(1).track_continuity, summary(2).track_continuity);
fprintf(fid, '| 丢轨次数 | %d | %d |\n', summary(1).lost_track_count, summary(2).lost_track_count);
fprintf(fid, '| 轨迹确认数 | %d | %d |\n', summary(1).confirmed_count, summary(2).confirmed_count);
fprintf(fid, '| 轨迹删除数 | %d | %d |\n\n', summary(1).deleted_count, summary(2).deleted_count);
fprintf(fid, '**结论：** %s\n\n', track_conclusion(summary));

fprintf(fid, '## 5 JPDA 分析\n\n');
fprintf(fid, '| 指标 | A (post) | B (L_softmax) | Δ |\n|------|----------|---------------|---|\n');
fprintf(fid, '| 平均门内量测数 | %.3f | %.3f | %+.3f |\n', ...
    summary(1).gate_meas_mean, summary(2).gate_meas_mean, ...
    summary(2).gate_meas_mean - summary(1).gate_meas_mean);
fprintf(fid, '| β₀ 平均值 | %.4f | %.4f | %+.4f |\n', ...
    summary(1).beta0_mean, summary(2).beta0_mean, summary(2).beta0_mean - summary(1).beta0_mean);
fprintf(fid, '| β_d 最大值均值 | %.4f | %.4f | %+.4f |\n', ...
    summary(1).beta_d_max_mean, summary(2).beta_d_max_mean, ...
    summary(2).beta_d_max_mean - summary(1).beta_d_max_mean);
fprintf(fid, '| 关联失败次数 | %d | %d | %d |\n', ...
    summary(1).assoc_fail_count, summary(2).assoc_fail_count, ...
    summary(2).assoc_fail_count - summary(1).assoc_fail_count);
fprintf(fid, '| 门控拒绝率 | %.4f | %.4f | %+.4f |\n\n', ...
    summary(1).gate_reject_rate, summary(2).gate_reject_rate, ...
    summary(2).gate_reject_rate - summary(1).gate_reject_rate);
fprintf(fid, '**结论：** %s\n\n', jpda_conclusion(summary));

fprintf(fid, '## 6 粒子分布分析\n\n');
fprintf(fid, '![粒子云对比](../figures/ab_particle_snapshots.png)\n\n');
fprintf(fid, '## 7 风险评估\n\n');
fprintf(fid, '| 问题 | 评估 |\n|------|------|\n');
fprintf(fid, '| 是否提高概率一致性？ | %s |\n', risk_prob_consistency(summary));
fprintf(fid, '| 是否改善粒子退化？ | %s |\n', risk_degeneracy(summary));
fprintf(fid, '| 是否降低跟踪精度？ | %s |\n', risk_accuracy(summary));
fprintf(fid, '| 是否影响 IMM 模型识别？ | %s |\n', risk_imm(summary));
fprintf(fid, '| 是否影响 JPDA 关联？ | %s |\n\n', risk_jpda(summary));

rec = decide_recommendation(summary);
fprintf(fid, '## 8 最终建议\n\n');
fprintf(fid, '**%s**\n\n', rec.text);
fprintf(fid, '%s\n', rec.rationale);

fclose(fid);
end


function t = ess_conclusion(s)
d_mean = s(2).ess_mean - s(1).ess_mean;
d_below = s(1).ess_below_pct - s(2).ess_below_pct;
if d_mean > 10 && d_below > 1
    t = 'B 组 ESS 显著更高，低 ESS 帧比例下降，粒子退化明显减轻。';
elseif d_mean > 0
    t = 'B 组 ESS 略有提升，粒子退化有所缓解。';
elseif d_mean < -10
    t = 'B 组 ESS 反而下降，需关注声学似然是否过尖。';
else
    t = '两组 ESS 差异不大。';
end
end


function t = imm_conclusion(s)
dm = abs(s(1).mu1_mean - 0.5) - abs(s(2).mu1_mean - 0.5);
if s(2).model_switch_count > s(1).model_switch_count * 1.5
    t = 'B 组模型切换增多，需关注机动识别稳定性。';
elseif dm > 0.05
    t = 'A 组更接近 μ≈0.5 僵死区；B 组模型概率分化更清晰。';
else
    t = '两组 IMM 行为接近，声学似然切换未显著改变模型识别。';
end
end


function t = track_conclusion(s)
if isnan(s(1).pos_rmse) || isnan(s(2).pos_rmse)
    t = '无有效 RMSE（轨迹未确认），请参考连续率与丢轨统计。';
elseif s(2).pos_rmse <= s(1).pos_rmse * 1.05
    t = 'B 组精度未劣于 A 组（RMSE 持平或更好）。';
elseif s(2).pos_rmse > s(1).pos_rmse * 1.15
    t = 'B 组 RMSE 明显劣于 A 组，存在精度代价。';
else
    t = '两组 RMSE 差异在 15% 以内，精度影响有限。';
end
end


function t = jpda_conclusion(s)
d = abs(s(2).gate_meas_mean - s(1).gate_meas_mean) + ...
    abs(s(2).beta0_mean - s(1).beta0_mean);
if d < 0.01 && s(1).assoc_fail_count == s(2).assoc_fail_count
    t = 'JPDA 层指标几乎不变，符合"关联仅依赖雷达量测"的理论预期。';
elseif d < 0.05
    t = 'JPDA 指标有微小数值差，来自轨迹状态反馈，关联层总体稳定。';
else
    t = 'JPDA 指标出现可观差异，需排查轨迹状态是否间接影响门控。';
end
end


function t = risk_prob_consistency(~)
t = '是。B 组仅使用当前帧 L_k，避免 post 历史信息与 PF 时间递推双重计入。';
end


function t = risk_degeneracy(s)
if s(2).ess_mean > s(1).ess_mean && s(2).total_resamples <= s(1).total_resamples
    t = '是。ESS 提升且重采样减少/持平。';
elseif s(2).ess_mean >= s(1).ess_mean
    t = '部分改善。ESS 提升但重采样未同步下降。';
else
    t = '未见改善，或略有退化。';
end
end


function t = risk_accuracy(s)
if isnan(s(2).pos_rmse), t = '无法判定（无 RMSE）。'; return; end
if s(2).pos_rmse <= s(1).pos_rmse * 1.05
    t = '否。B 组精度未下降。';
elseif s(2).pos_rmse > s(1).pos_rmse * 1.1
    t = '是。B 组 RMSE 上升超过 10%。';
else
    t = '轻微。RMSE 变化在 10% 以内。';
end
end


function t = risk_imm(s)
if abs(s(2).mu1_mean - 0.5) < abs(s(1).mu1_mean - 0.5)
    t = '否/略有改善。B 组未加剧 μ≈0.5 僵死。';
else
    t = '可能有影响，需结合模型切换次数判断。';
end
end


function t = risk_jpda(s)
d = max(abs(s(2).gate_meas_mean - s(1).gate_meas_mean), ...
          abs(s(2).beta0_mean - s(1).beta0_mean));
if d < 0.02, t = '否。JPDA 统计量基本一致。'; else, t = '有间接影响，但通常来自轨迹状态差异。'; end
end


function rec = decide_recommendation(s)
% 只能 A / B / C 三选一
score = 0;
if s(2).ess_mean > s(1).ess_mean + 5, score = score + 2; end
if s(2).ess_below_pct < s(1).ess_below_pct - 1, score = score + 1; end
if s(2).total_resamples <= s(1).total_resamples, score = score + 1; end
if ~isnan(s(2).pos_rmse) && s(2).pos_rmse <= s(1).pos_rmse * 1.1, score = score + 2; end
if abs(s(2).gate_meas_mean - s(1).gate_meas_mean) < 0.05, score = score + 1; end

if score >= 5
    rec.text = 'A. 建议保留修改（L_softmax）';
    rec.rationale = 'B 组在 ESS/退化与精度上均不劣于 A，且概率一致性更符合贝叶斯滤波框架。';
elseif score <= 1 && s(2).pos_rmse > s(1).pos_rmse * 1.15
    rec.text = 'B. 建议回滚修改（恢复 post）';
    rec.rationale = 'B 组未带来 ESS 改善，且跟踪精度显著下降。';
else
    rec.text = 'C. 需要进一步实验';
    rec.rationale = '两组差异不够显著或存在精度/退化权衡，建议扩展场景或调参前再验证。';
end
end


function s = bool_str(v)
if v, s = '是'; else, s = '否'; end
end


function out = iif(cond, a, b)
if cond, out = a; else, out = b; end
end
