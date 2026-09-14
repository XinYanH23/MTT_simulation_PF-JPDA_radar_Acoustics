%% 实验一补充：声学知情删除 + 场辅助确认（方法 A / B 对照）
% 默认 config 开关为关；本脚本临时打开，不覆盖 exp1 主结果。
clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'IMM'));
addpath(fullfile(this_dir, 'JPDA'));
addpath(fullfile(this_dir, 'EKF'));
addpath(fullfile(this_dir, 'PF'));
addpath(fullfile(this_dir, 'TrackManagement'));
addpath(fullfile(this_dir, 'Utilities'));
addpath(fullfile(this_dir, 'Evaluation'));
addpath(fullfile(this_dir, 'Scenario'));
addpath(fullfile(this_dir, 'matlab'));

out_dir = fullfile(this_dir, 'experiments', 'exp1_ac_lifecycle');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp1_life');
env_out = getenv('TRACKER_OUT_DIR');
env_tex = getenv('TRACKER_TEX_DIR');
if ~isempty(env_out); out_dir = env_out; end
if ~isempty(env_tex); tex_dir = env_tex; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

MC_DEBUG = false;
n_mc     = 25;
if strcmp(getenv('TRACKER_MC_DEBUG'), '1'); MC_DEBUG = true; end
env_n = str2double(getenv('TRACKER_N_MC'));
if ~isnan(env_n) && env_n >= 1; n_mc = env_n; end
[seeds, n_mc, is_debug] = mc_resolve_seeds(MC_DEBUG, n_mc); %#ok<ASGLU>

fprintf('===== 声学生命周期：知情删除 + 场辅助确认 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time);

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));

cases = { ...
    struct('key','A_life', 'mode','legacy_likelihood', 'label','方法A + 知情生命周期'), ...
    struct('key','B_life', 'mode','jpda_posterior',    'label','方法B + 知情生命周期') };

opts = struct('enable_acoustic', true, 'verbose', false, 'method', 'softmax');
packs = struct();
cache_mat = fullfile(out_dir, 'exp1_ac_lifecycle.mat');
if exist(cache_mat, 'file')
    old = load(cache_mat, 'packs');
    if isfield(old, 'packs'); packs = old.packs; end
end

for ci = 1:numel(cases)
    cs = cases{ci};
    if isfield(packs, cs.key) && isfield(packs.(cs.key), 'agg')
        fprintf('\n########## %s  跳过（缓存） ##########\n', cs.label);
        continue
    end
    fprintf('\n########## %s  (%d seeds) ##########\n', cs.label, n_mc);
    cfg = config_tracker();
    cfg.use_pf = true;
    cfg.save_mat = false;
    cfg.jpda_radar_only = true;
    cfg.birth_radar_only = true;
    cfg.pf_update_mode = cs.mode;
    cfg.ac_informed_delete = true;
    cfg.ac_aided_confirm   = true;
    packs.(cs.key) = mc_run_case(allDetections, truth, time, cfg, ac_cfg, opts, seeds);
    packs.(cs.key).label = cs.label;
    packs.(cs.key).mode  = cs.mode;
    save(cache_mat, 'packs', 'cases', 'truth', 'time', 'seeds', 'n_mc', 'is_debug', '-v7.3');
end

% 基线：已完成的实验一 Field（只读）
base = struct();
base.A = try_load_field(fullfile(this_dir, 'methods_v1_legacy', ...
    'experiments', 'exp1_acoustic_ablation', 'exp1_results.mat'), ...
    fullfile(this_dir, 'experiments', 'exp1_acoustic_ablation', 'exp1_results.mat'));
base.B = try_load_field(fullfile(this_dir, 'methods_v2_bayes', ...
    'experiments', 'exp1_acoustic_ablation', 'exp1_results.mat'), '');

fid = fopen(fullfile(out_dir, 'exp1_ac_lifecycle.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验一补充：声学知情删除 + 场辅助确认\n\n');
fprintf(fid, '删除：场弱且连续漏检≥8 → 提前删；场强则保活至 30 帧。\n');
fprintf(fid, '确认：场强 2/5，场弱仍 3/5。诞生仍仅未关联雷达点 + η_a。\n');
fprintf(fid, '种子 %s。\n\n', mat2str(seeds));
fprintf(fid, '| 方案 | 峰值确认 | 时均确认 | OSPA | 虚警 | 基数误差 | 中位寿命 | 碎裂 |\n');
fprintf(fid, '|------|----------|----------|------|------|----------|----------|------|\n');
rows = {};
if isstruct(base.A)
    rows(end+1,:) = {'方法A 基线 Field', base.A}; %#ok<AGROW>
end
rows(end+1,:) = {'方法A + 生命周期', packs.A_life}; %#ok<AGROW>
if isstruct(base.B)
    rows(end+1,:) = {'方法B 基线 Field', base.B}; %#ok<AGROW>
end
rows(end+1,:) = {'方法B + 生命周期', packs.B_life}; %#ok<AGROW>
for i = 1:size(rows,1)
    a = rows{i,2}.agg;
    mn = mean_n_est_stats(rows{i,2});
    fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s | %s |\n', rows{i,1}, ...
        mc_pm_md(a.peak_confirmed,'f1'), mc_pm_md(mn,'f2'), ...
        mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
        mc_pm_md(a.card_err_mean), mc_pm_md(a.median_track_life,'f1'), ...
        mc_pm_md(a.fragmentation_count,'f1'));
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp1_ac_lifecycle.md'), fullfile(tex_dir, 'exp1_ac_lifecycle.md'));

t_sec = time(:)';
if max(t_sec) < 1; t_sec = (0:N-1) * 0.1; end
n_true = packs.A_life.agg.n_true_k;
fig = figure('Color','w','Position',[60 60 920 420]);
plot(t_sec, n_true, 'k-', 'LineWidth', 1.8); hold on;
if isstruct(base.A)
    mc_plot_series_ci(t_sec, base.A.agg.n_est_k, [0.75 0.55 0.20], 25);
end
mc_plot_series_ci(t_sec, packs.A_life.agg.n_est_k, [0.80 0.25 0.20], 25);
if isstruct(base.B)
    mc_plot_series_ci(t_sec, base.B.agg.n_est_k, [0.45 0.70 0.85], 25);
end
mc_plot_series_ci(t_sec, packs.B_life.agg.n_est_k, [0.15 0.35 0.70], 25);
grid on; box on;
xlabel('Time (s)'); ylabel('Number of confirmed tracks');
leg = {'N_{true}'};
if isstruct(base.A); leg{end+1} = 'A baseline'; end %#ok<AGROW>
leg{end+1} = 'A + lifecycle';
if isstruct(base.B); leg{end+1} = 'B baseline'; end %#ok<AGROW>
leg{end+1} = 'B + lifecycle';
legend(leg, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
exportgraphics(fig, fullfile(out_dir, 'fig_exp1_ac_lifecycle_N.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(tex_dir, 'fig_exp1_ac_lifecycle_N.png'), 'Resolution', 300);
fprintf('\n===== 声学生命周期扫描完成 =====\n');

function pack = try_load_field(p1, p2)
pack = [];
for p = {p1, p2}
    if isempty(p{1}) || ~exist(p{1}, 'file'); continue; end
    S = load(p{1});
    if isfield(S, 'mc_B')
        pack = S.mc_B;
        return
    end
end
end

function s = mean_n_est_stats(pack)
x = nan(numel(pack.per_seed), 1);
for i = 1:numel(pack.per_seed)
    x(i) = mean(pack.per_seed{i}.n_est_k);
end
s = mc_scalar_stats(x);
end
