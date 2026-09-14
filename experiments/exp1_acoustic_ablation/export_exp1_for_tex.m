%% EXPORT_EXP1_FOR_TEX  从 exp1_results.mat（蒙特卡洛）重绘并写表
clear; clc; close all;

exp_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(exp_dir));
addpath(fullfile(root, 'Evaluation'));
tex_dir = fullfile(root, 'report', 'figures', 'ch3_exp1');
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

S = load(fullfile(exp_dir, 'exp1_results.mat'));
assert(isfield(S, 'mc_A') && isfield(S, 'mc_B'), ...
    'exp1_results.mat 需含 mc_A/mc_B（请先跑 run_exp1_acoustic_ablation.m）');
A = S.mc_A.agg; B = S.mc_B.agg;
time = S.time; cfg = S.cfg;
n_mc = S.n_mc;
is_debug = false;
if isfield(S, 'is_debug'); is_debug = S.is_debug; end
t_sec = time(:)';
if max(t_sec) < 1; t_sec = (0:numel(time)-1)' * cfg.dt; end

colA = [0.75 0.25 0.25];
colB = [0.15 0.35 0.70];

fig1 = figure('Color','w','Position',[80 80 900 380]);
mc_plot_series_ci(t_sec, A.ospa_k, colA, 25); hold on;
mc_plot_series_ci(t_sec, B.ospa_k, colB, 25);
grid on; box on;
xlabel('Time (s)'); ylabel('OSPA (m)');
legend({sprintf('Radar-only (%s)', mc_pm_md(A.ospa_mean)), ...
        sprintf('Radar+Acoustic (%s)', mc_pm_md(B.ospa_mean))}, ...
        'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig1, tex_dir, 'fig_exp1_ospa');

fig2 = figure('Color','w','Position',[80 80 900 380]);
plot(t_sec, A.n_true_k, 'k-', 'LineWidth', 1.6); hold on;
mc_plot_series_ci(t_sec, A.n_est_k, colA, 25);
mc_plot_series_ci(t_sec, B.n_est_k, colB, 25);
grid on; box on;
xlabel('Time (s)'); ylabel('Number of tracks');
legend({'$N_{\mathrm{true}}$', 'Radar-only $\hat{N}$', 'Radar+Acoustic $\hat{N}$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig2, tex_dir, 'fig_exp1_cardinality');

fig3 = figure('Color','w','Position',[80 80 480 380]);
sA = scale_pct(A.false_confirm_ratio);
sB = scale_pct(B.false_confirm_ratio);
mc_bar_with_ci({sA, sB}, [colA; colB], {'Radar-only', 'Radar+Acoustic'});
ylabel('False confirmed track ratio (\%)');
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig3, tex_dir, 'fig_exp1_false_confirm');

fig4 = figure('Color','w','Position',[80 80 900 380]);
subplot(1,2,1);
mc_bar_with_ci({A.mean_track_life, B.mean_track_life}, [colA; colB], ...
    {'Radar-only', 'Radar+Acoustic'});
ylabel('Mean track life (frames)'); grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
subplot(1,2,2);
mc_bar_with_ci({A.median_track_life, B.median_track_life}, [colA; colB], ...
    {'Radar-only', 'Radar+Acoustic'});
ylabel('Median track life (frames)'); grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig4, tex_dir, 'fig_exp1_track_life');

% table
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
fid = fopen(fullfile(tex_dir, 'tab_exp1_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp1 MC table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{实验一：Radar-only 与 Radar+Acoustic Field 对比（$\\gamma=0.3$，%s）}\n', cap_note);
fprintf(fid, '  \\label{tab:exp1_metrics}\n');
fprintf(fid, '  \\begin{tabular}{lcc}\n    \\toprule\n');
fprintf(fid, '    指标 & Radar-only & Radar+Acoustic \\\\\n    \\midrule\n');
fprintf(fid, '    OSPA 均值 (m) & %s & %s \\\\\n', mc_pm_tex(A.ospa_mean), mc_pm_tex(B.ospa_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{loc}}$ (m) & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_loc_mean), mc_pm_tex(B.ospa_loc_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{card}}$ (m) & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_card_mean), mc_pm_tex(B.ospa_card_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{assign}}$ (m) & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_assign_mean), mc_pm_tex(B.ospa_assign_mean));
fprintf(fid, '    OSPA RMSE (m) & %s & %s \\\\\n', mc_pm_tex(A.ospa_rmse), mc_pm_tex(B.ospa_rmse));
fprintf(fid, '    基数误差均值 & %s & %s \\\\\n', mc_pm_tex(A.card_err_mean), mc_pm_tex(B.card_err_mean));
fprintf(fid, '    峰值确认航迹数 & %s & %s \\\\\n', mc_pm_tex(A.peak_confirmed,'f1'), mc_pm_tex(B.peak_confirmed,'f1'));
fprintf(fid, '    虚警确认比例 & %s & %s \\\\\n', mc_pm_tex(A.false_confirm_ratio,'pct'), mc_pm_tex(B.false_confirm_ratio,'pct'));
fprintf(fid, '    中位航迹寿命 (帧) & %s & %s \\\\\n', mc_pm_tex(A.median_track_life,'f1'), mc_pm_tex(B.median_track_life,'f1'));
fprintf(fid, '    碎裂次数 & %s & %s \\\\\n', mc_pm_tex(A.fragmentation_count,'f1'), mc_pm_tex(B.fragmentation_count,'f1'));
fprintf(fid, '    ID Switch & %s & %s \\\\\n', mc_pm_tex(A.id_switches,'f1'), mc_pm_tex(B.id_switches,'f1'));
fprintf(fid, '    位置 RMSE (m，附带) & %s & %s \\\\\n', mc_pm_tex(A.pos_rmse), mc_pm_tex(B.pos_rmse));
fprintf(fid, '    \\bottomrule\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);

fprintf('TeX figures/table exported to: %s\n', tex_dir);

function s = scale_pct(s)
s.mean = 100*s.mean; s.std = 100*s.std;
s.ci_lo = 100*s.ci_lo; s.ci_hi = 100*s.ci_hi;
end

function export_both(fig, out_dir, name)
png = fullfile(out_dir, [name '.png']);
pdf = fullfile(out_dir, [name '.pdf']);
exportgraphics(fig, png, 'Resolution', 300);
exportgraphics(fig, pdf, 'ContentType', 'vector');
fprintf('  saved %s .{png,pdf}\n', name);
end
