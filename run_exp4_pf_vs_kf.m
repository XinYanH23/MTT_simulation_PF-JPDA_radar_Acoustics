%% =========================================================================
%  run_exp4_pf_vs_kf.m
%  实验四：同一声学场 + 同一乘性更新下 IMM-KF vs IMM-PF（25 组独立种子；rng(42) 仅调试）
%% =========================================================================
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

out_dir = fullfile(this_dir, 'experiments', 'exp4_pf_vs_kf');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp4');
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

fprintf('===== 实验四：加载场景与探测 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time);

cfg0 = config_tracker();
cfg0.save_mat = false;
cfg0.jpda_radar_only = true;
cfg0.birth_radar_only = true;
cfg0.pf_update_mode = 'legacy_likelihood';
env_mode = getenv('TRACKER_PF_UPDATE_MODE');
if ~isempty(env_mode); cfg0.pf_update_mode = env_mode; end
fprintf('  pf_update_mode=%s\n  out_dir=%s\n', cfg0.pf_update_mode, out_dir);
cfg0.pf_ac_gamma_mode = 'fixed';
cfg0.pf_ac_gamma_fixed = 0.3;

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));

cases = { ...
    struct('key','KF_Field',  'use_pf', false, 'mode','field', 'label','IMM-KF + Field'), ...
    struct('key','PF_Field',  'use_pf', true,  'mode','field', 'label','IMM-PF + Field'), ...
    struct('key','PF_Radar',  'use_pf', true,  'mode','none',  'label','IMM-PF Radar-only') };

packs = struct();
for ci = 1:numel(cases)
    cs = cases{ci};
    fprintf('\n########## %s  (%d seeds) ##########\n', cs.label, n_mc);
    cfg = cfg0;
    cfg.use_pf = cs.use_pf;
    opts = struct('acoustic_mode', cs.mode, 'verbose', false, 'method', 'softmax');
    packs.(cs.key) = mc_run_case(allDetections, truth, time, cfg, ac_cfg, opts, seeds);
end

save(fullfile(out_dir, 'exp4_results.mat'), ...
    'packs', 'cases', 'truth', 'time', 'cfg0', 'seeds', 'n_mc', 'is_debug', '-v7.3');
mc_write_seed_csv(fullfile(out_dir, 'exp4_per_seed.csv'), ...
    {packs.KF_Field, packs.PF_Field, packs.PF_Radar}, ...
    {'KF_Field', 'PF_Field', 'PF_Radar'});

keys = {'KF_Field','PF_Field','PF_Radar'};
t_sec = time(:); if max(t_sec) < 1; t_sec = (0:N-1)' * cfg0.dt; end

fid = fopen(fullfile(out_dir, 'exp4_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验四：IMM-KF vs IMM-PF（蒙特卡洛）\n\n');
fprintf(fid, '种子：主实验 %d 组 %s；rng(42) 仅为调试。\n\n', n_mc, mat2str(seeds));
fprintf(fid, '| 方案 | OSPA | 虚警确认 | 基数误差 | 中位寿命 | 碎裂 | IDS |\n');
fprintf(fid, '|------|------|----------|----------|----------|------|-----|\n');
for i = 1:numel(keys)
    a = packs.(keys{i}).agg;
    fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s |\n', cases{i}.label, ...
        mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
        mc_pm_md(a.card_err_mean), mc_pm_md(a.median_track_life,'f1'), ...
        mc_pm_md(a.fragmentation_count,'f1'), mc_pm_md(a.id_switches,'f1'));
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp4_summary.md'), fullfile(tex_dir, 'exp4_summary.md'));

col = [0.55 0.35 0.15; 0.15 0.35 0.70; 0.75 0.25 0.25];
fig1 = figure('Color','w','Position',[80 80 900 380]);
hold on;
for i = 1:numel(keys)
    mc_plot_series_ci(t_sec, packs.(keys{i}).agg.ospa_k, col(i,:), 25);
end
grid on; box on;
xlabel('Time (s)'); ylabel('OSPA (m)');
legend({cases{1}.label, cases{2}.label, cases{3}.label}, 'Location', 'northeast', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig1, out_dir, tex_dir, 'fig_exp4_ospa');

fig2 = figure('Color','w','Position',[80 80 560 380]);
slist = {scale_pct(packs.KF_Field.agg.false_confirm_ratio), ...
         scale_pct(packs.PF_Field.agg.false_confirm_ratio), ...
         scale_pct(packs.PF_Radar.agg.false_confirm_ratio)};
mc_bar_with_ci(slist, col, {'KF+Field', 'PF+Field', 'PF Radar'});
ylabel('False confirmed track ratio (%)');
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig2, out_dir, tex_dir, 'fig_exp4_false_confirm');

fig3 = figure('Color','w','Position',[80 80 720 380]);
plot(t_sec, packs.KF_Field.agg.n_true_k, 'k-', 'LineWidth', 1.5); hold on;
mc_plot_series_ci(t_sec, packs.KF_Field.agg.n_est_k, col(1,:), 25);
mc_plot_series_ci(t_sec, packs.PF_Field.agg.n_est_k, col(2,:), 25);
grid on; box on;
xlabel('Time (s)'); ylabel('Number of tracks');
legend({'N_{true}', 'KF+Field', 'PF+Field'}, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig3, out_dir, tex_dir, 'fig_exp4_cardinality');

write_exp4_tex_table(tex_dir, packs, n_mc, is_debug);
fprintf('\n===== 实验四完成（n_mc=%d） =====\n', n_mc);

function write_exp4_tex_table(tex_dir, packs, n_mc, is_debug)
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
K = packs.KF_Field.agg; P = packs.PF_Field.agg; R = packs.PF_Radar.agg;
fid = fopen(fullfile(tex_dir, 'tab_exp4_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp4 MC table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{实验四：对等场输入与对等乘性更新下的 KF/PF 对比（$\\gamma=0.3$，%s）}\n', cap_note);
fprintf(fid, '  \\label{tab:exp4_metrics}\n');
fprintf(fid, '  \\begin{tabular}{lccc}\n    \\toprule\n');
fprintf(fid, '    指标 & KF+Field & PF+Field & PF Radar-only \\\\\n    \\midrule\n');
fprintf(fid, '    OSPA 均值 (m) & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.ospa_mean), mc_pm_tex(P.ospa_mean), mc_pm_tex(R.ospa_mean));
fprintf(fid, '    虚警确认比例 & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.false_confirm_ratio,'pct'), mc_pm_tex(P.false_confirm_ratio,'pct'), mc_pm_tex(R.false_confirm_ratio,'pct'));
fprintf(fid, '    基数误差均值 & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.card_err_mean), mc_pm_tex(P.card_err_mean), mc_pm_tex(R.card_err_mean));
fprintf(fid, '    中位航迹寿命 (帧) & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.median_track_life,'f1'), mc_pm_tex(P.median_track_life,'f1'), mc_pm_tex(R.median_track_life,'f1'));
fprintf(fid, '    碎裂次数 & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.fragmentation_count,'f1'), mc_pm_tex(P.fragmentation_count,'f1'), mc_pm_tex(R.fragmentation_count,'f1'));
fprintf(fid, '    ID Switch & %s & %s & %s \\\\\n', ...
    mc_pm_tex(K.id_switches,'f1'), mc_pm_tex(P.id_switches,'f1'), mc_pm_tex(R.id_switches,'f1'));
fprintf(fid, '    \\bottomrule\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);
end

function s = scale_pct(s)
s.mean = 100*s.mean; s.std = 100*s.std;
s.ci_lo = 100*s.ci_lo; s.ci_hi = 100*s.ci_hi;
end

function export_both(fig, out_dir, tex_dir, name)
exportgraphics(fig, fullfile(out_dir, [name '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, [name '.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(tex_dir, [name '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(tex_dir, [name '.pdf']), 'ContentType', 'vector');
fprintf('  saved %s\n', name);
end
