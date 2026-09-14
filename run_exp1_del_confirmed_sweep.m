%% 实验一补充：Confirmed→Deleted 连续漏检门限扫描
% 不改 config_tracker 默认 30，不覆盖 exp1 主结果。
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

out_dir = fullfile(this_dir, 'experiments', 'exp1_del_confirmed_sweep');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp1_del');
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

del_list = [15, 20, 25];
env_d = getenv('TRACKER_DEL_LIST');
if ~isempty(env_d)
    del_list = str2num(env_d); %#ok<ST2NM>
    del_list = unique(del_list(:)');
end

fprintf('===== 实验一补充：del_confirmed 扫描 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time);

cfg0 = config_tracker();
cfg0.use_pf = true;
cfg0.save_mat = false;
cfg0.jpda_radar_only = true;
cfg0.birth_radar_only = true;
cfg0.pf_update_mode = 'legacy_likelihood';
env_mode = getenv('TRACKER_PF_UPDATE_MODE');
if ~isempty(env_mode); cfg0.pf_update_mode = env_mode; end
fprintf('  pf_update_mode=%s  del_list=%s  n_mc=%d\n', ...
    cfg0.pf_update_mode, mat2str(del_list), n_mc);
fprintf('  out_dir=%s\n', out_dir);

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));

modes = { ...
    struct('key','Radar', 'enable', false, 'eta', 0,    'label','Radar-only'), ...
    struct('key','Field', 'enable', true,  'eta', [],   'label','Radar+Field') };

packs = struct();
cache_mat = fullfile(out_dir, 'exp1_del_results.mat');
if exist(cache_mat, 'file')
    old = load(cache_mat, 'packs');
    if isfield(old, 'packs'); packs = old.packs; end
    fprintf('  已加载缓存\n');
end

for di = 1:numel(del_list)
    dlt = del_list(di);
    for mi = 1:numel(modes)
        md = modes{mi};
        key = sprintf('%s_d%d', md.key, dlt);
        if isfield(packs, key) && isfield(packs.(key), 'agg')
            fprintf('\n########## %s  del=%d  跳过（缓存） ##########\n', md.label, dlt);
            continue
        end
        fprintf('\n########## %s  del_confirmed=%d  (%d seeds) ##########\n', ...
            md.label, dlt, n_mc);
        cfg = cfg0;
        cfg.del_confirmed = dlt;
        if ~isempty(md.eta); cfg.birth_acoustic_eta = md.eta; end
        opts = struct('enable_acoustic', md.enable, 'verbose', false, 'method', 'softmax');
        packs.(key) = mc_run_case(allDetections, truth, time, cfg, ac_cfg, opts, seeds);
        packs.(key).del_confirmed = dlt;
        packs.(key).label = md.label;
        save(cache_mat, 'packs', 'del_list', 'truth', 'time', 'cfg0', ...
            'seeds', 'n_mc', 'is_debug', '-v7.3');
    end
end

% 并入同方法已完成的实验一 del=30（只读对照）。勿跨方法混用。
ref30 = getenv('TRACKER_REF30_MAT');
if isempty(ref30)
    ref30 = fullfile(fileparts(out_dir), 'exp1_acoustic_ablation', 'exp1_results.mat');
end
has30 = false;
if exist(ref30, 'file')
    R = load(ref30, 'mc_A', 'mc_B');
    if isfield(R, 'mc_A') && isfield(R, 'mc_B')
        packs.Radar_d30 = R.mc_A;
        packs.Radar_d30.del_confirmed = 30;
        packs.Radar_d30.label = 'Radar-only';
        packs.Field_d30 = R.mc_B;
        packs.Field_d30.del_confirmed = 30;
        packs.Field_d30.label = 'Radar+Field';
        has30 = true;
        fprintf('  已并入实验一 del=30 对照：%s\n', ref30);
    end
end

del_plot = del_list;
if has30; del_plot = unique([del_list, 30]); end

save(cache_mat, 'packs', 'del_list', 'del_plot', 'truth', 'time', 'cfg0', ...
    'seeds', 'n_mc', 'is_debug', 'has30', '-v7.3');

fid = fopen(fullfile(out_dir, 'exp1_del_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验一补充：del_confirmed 扫描\n\n');
fprintf(fid, 'pf_update_mode=`%s`；种子 %s；默认 30 未改，本表为对照扫描。\n\n', ...
    cfg0.pf_update_mode, mat2str(seeds));
fprintf(fid, '| del | 方案 | 峰值确认航迹 | 时均确认航迹 | OSPA | 虚警确认 | 基数误差 | 中位寿命 | 碎裂 |\n');
fprintf(fid, '|----:|------|--------------|--------------|------|----------|----------|----------|------|\n');
for dlt = del_plot
    for mi = 1:numel(modes)
        key = sprintf('%s_d%d', modes{mi}.key, dlt);
        if ~isfield(packs, key); continue; end
        a = packs.(key).agg;
        mn = mean_n_est_stats(packs.(key));
        fprintf(fid, '| %d | %s | %s | %s | %s | %s | %s | %s | %s |\n', dlt, ...
            modes{mi}.label, mc_pm_md(a.peak_confirmed,'f1'), mc_pm_md(mn,'f2'), ...
            mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
            mc_pm_md(a.card_err_mean), mc_pm_md(a.median_track_life,'f1'), ...
            mc_pm_md(a.fragmentation_count,'f1'));
    end
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp1_del_summary.md'), fullfile(tex_dir, 'exp1_del_summary.md'));

colR = [0.75 0.25 0.25];
colF = [0.15 0.35 0.70];
sR_ospa = cell(1, numel(del_plot)); sF_ospa = sR_ospa;
sR_fa = sR_ospa; sF_fa = sR_ospa;
sR_peak = sR_ospa; sF_peak = sR_ospa;
sR_nbar = sR_ospa; sF_nbar = sR_ospa;
for i = 1:numel(del_plot)
    sR_ospa{i} = packs.(sprintf('Radar_d%d', del_plot(i))).agg.ospa_mean;
    sF_ospa{i} = packs.(sprintf('Field_d%d', del_plot(i))).agg.ospa_mean;
    sR_fa{i}   = scale_pct(packs.(sprintf('Radar_d%d', del_plot(i))).agg.false_confirm_ratio);
    sF_fa{i}   = scale_pct(packs.(sprintf('Field_d%d', del_plot(i))).agg.false_confirm_ratio);
    sR_peak{i} = packs.(sprintf('Radar_d%d', del_plot(i))).agg.peak_confirmed;
    sF_peak{i} = packs.(sprintf('Field_d%d', del_plot(i))).agg.peak_confirmed;
    sR_nbar{i} = mean_n_est_stats(packs.(sprintf('Radar_d%d', del_plot(i))));
    sF_nbar{i} = mean_n_est_stats(packs.(sprintf('Field_d%d', del_plot(i))));
end

fig = figure('Color','w','Position',[60 60 980 720]);
subplot(2,2,1);
mc_errorbar_xy(del_plot, sR_peak, colR, 's'); hold on;
mc_errorbar_xy(del_plot, sF_peak, colF, 'o');
yline(5, 'k--', 'N_{true} peak', 'LabelHorizontalAlignment', 'left');
grid on; box on; xlabel('del_{confirmed} (frames)'); ylabel('Peak confirmed tracks');
legend({'Radar-only','Radar+Field'}, 'Location','best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', del_plot);
subplot(2,2,2);
mc_errorbar_xy(del_plot, sR_nbar, colR, 's'); hold on;
mc_errorbar_xy(del_plot, sF_nbar, colF, 'o');
grid on; box on; xlabel('del_{confirmed} (frames)'); ylabel('Time-mean confirmed tracks');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', del_plot);
subplot(2,2,3);
mc_errorbar_xy(del_plot, sR_ospa, colR, 's'); hold on;
mc_errorbar_xy(del_plot, sF_ospa, colF, 'o');
grid on; box on; xlabel('del_{confirmed} (frames)'); ylabel('OSPA mean (m)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', del_plot);
subplot(2,2,4);
mc_errorbar_xy(del_plot, sR_fa, colR, 's'); hold on;
mc_errorbar_xy(del_plot, sF_fa, colF, 'o');
grid on; box on; xlabel('del_{confirmed} (frames)'); ylabel('False confirmed (%)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', del_plot);
exportgraphics(fig, fullfile(out_dir, 'fig_exp1_del_sweep.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, 'fig_exp1_del_sweep.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(tex_dir, 'fig_exp1_del_sweep.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(tex_dir, 'fig_exp1_del_sweep.pdf'), 'ContentType', 'vector');

t_sec = time(:)';
if max(t_sec) < 1; t_sec = (0:N-1) * cfg0.dt; end
n_true = packs.(sprintf('Field_d%d', del_plot(1))).agg.n_true_k;
colsF = [0.70 0.45 0.10; 0.20 0.55 0.30; 0.15 0.35 0.70; 0.45 0.20 0.55];
figc = figure('Color','w','Position',[60 60 920 420]);
plot(t_sec, n_true, 'k-', 'LineWidth', 1.8); hold on;
leg = {'N_{true}'};
for i = 1:numel(del_plot)
    key = sprintf('Field_d%d', del_plot(i));
    if ~isfield(packs, key); continue; end
    ci = 1 + mod(i-1, size(colsF,1));
    mc_plot_series_ci(t_sec, packs.(key).agg.n_est_k, colsF(ci,:), 25);
    leg{end+1} = sprintf('Field del=%d', del_plot(i)); %#ok<AGROW>
end
grid on; box on;
xlabel('Time (s)'); ylabel('Number of confirmed tracks');
legend(leg, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
exportgraphics(figc, fullfile(out_dir, 'fig_exp1_del_cardinality.png'), 'Resolution', 300);
exportgraphics(figc, fullfile(out_dir, 'fig_exp1_del_cardinality.pdf'), 'ContentType', 'vector');
exportgraphics(figc, fullfile(tex_dir, 'fig_exp1_del_cardinality.png'), 'Resolution', 300);
exportgraphics(figc, fullfile(tex_dir, 'fig_exp1_del_cardinality.pdf'), 'ContentType', 'vector');

fprintf('\n===== del_confirmed 扫描完成 =====\n');

function s = scale_pct(s)
s.mean = 100*s.mean; s.std = 100*s.std;
s.ci_lo = 100*s.ci_lo; s.ci_hi = 100*s.ci_hi;
end

function s = mean_n_est_stats(pack)
x = nan(numel(pack.per_seed), 1);
for i = 1:numel(pack.per_seed)
    x(i) = mean(pack.per_seed{i}.n_est_k);
end
s = mc_scalar_stats(x);
end
