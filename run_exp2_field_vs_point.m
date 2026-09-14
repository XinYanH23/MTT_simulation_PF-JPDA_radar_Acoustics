%% =========================================================================
%  run_exp2_field_vs_point.m
%  实验二：声学概率场 vs 点量测对比
%
%  主实验：25 组独立种子；rng(42) 仅为调试。
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

out_dir = fullfile(this_dir, 'experiments', 'exp2_field_vs_point');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp2');
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

fprintf('===== 实验二：加载场景与探测 =====\n');
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
fprintf('  pf_update_mode=%s\n  out_dir=%s\n', cfg0.pf_update_mode, out_dir);
cfg0.pf_ac_gamma_mode = 'fixed';
cfg0.pf_ac_gamma_fixed = 0.3;

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));
fprintf('  帧=%d UAV=%d 声学站=%d γ=%.1f\n', N, numel(truth), size(data1.acousticPos,1), cfg0.pf_ac_gamma_fixed);

modes = { ...
    struct('name','Baseline-R',    'key','Baseline_R',    'acoustic_mode','none',        'label','Radar-only'), ...
    struct('name','Point-A',       'key','Point_A',       'acoustic_mode','point_map',   'label','Point (MAP)'), ...
    struct('name','Point-A-Multi', 'key','Point_A_Multi', 'acoustic_mode','point_multi', 'label','Point (Multi)'), ...
    struct('name','Field-A',       'key','Field_A',       'acoustic_mode','field',       'label','Field (Ours)') };

packs = struct();
for mi = 1:numel(modes)
    md = modes{mi};
    fprintf('\n########## %s [%s]  (%d seeds) ##########\n', md.name, md.acoustic_mode, n_mc);
    opts = struct('acoustic_mode', md.acoustic_mode, 'verbose', false, ...
                  'method', 'softmax', 'point_K', 3);
    packs.(md.key) = mc_run_case(allDetections, truth, time, cfg0, ac_cfg, opts, seeds);
end

save(fullfile(out_dir, 'exp2_results.mat'), ...
    'packs', 'modes', 'truth', 'time', 'cfg0', 'seeds', 'n_mc', 'is_debug', '-v7.3');
mc_write_seed_csv(fullfile(out_dir, 'exp2_per_seed.csv'), ...
    {packs.Baseline_R, packs.Point_A, packs.Point_A_Multi, packs.Field_A}, ...
    {'Baseline-R', 'Point-A', 'Point-A-Multi', 'Field-A'});

keys  = {'Baseline_R','Point_A','Point_A_Multi','Field_A'};
names = {'Baseline-R','Point-A','Point-A-Multi','Field-A'};
t_sec = time(:)';
if max(t_sec) < 1; t_sec = (0:N-1)' * cfg0.dt; end

write_exp2_summary(out_dir, tex_dir, packs, keys, names, cfg0, data1, seeds, n_mc);
plot_exp2_figures(out_dir, tex_dir, packs, keys, names, t_sec);
write_exp2_tex_table(tex_dir, packs, keys, n_mc, is_debug);

fprintf('\n===== 实验二完成（n_mc=%d） =====\n', n_mc);

function write_exp2_summary(out_dir, tex_dir, packs, keys, names, cfg, data1, seeds, n_mc)
fid = fopen(fullfile(out_dir, 'exp2_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验二：声学概率场 vs 点量测（蒙特卡洛）\n\n');
fprintf(fid, '| 项目 | 取值 |\n|------|------|\n');
fprintf(fid, '| 随机种子 | 主实验 %d 组 %s；rng(42) 仅为调试 |\n', n_mc, mat2str(seeds));
fprintf(fid, '| PF | N=%d，γ=%.1f |\n', cfg.pf_N, cfg.pf_ac_gamma_fixed);
fprintf(fid, '| 声学站 | %d |\n\n', size(data1.acousticPos,1));
fprintf(fid, '## 结果（均值 ± 标准差）\n\n| 指标 |');
for i=1:4, fprintf(fid, ' %s |', names{i}); end
fprintf(fid, '\n|------|');
for i=1:4, fprintf(fid, '------|'); end
fprintf(fid, '\n');
rows = { ...
    'OSPA 均值 (m)',     @(a) mc_pm_md(a.ospa_mean); ...
    'OSPA loc (m)',      @(a) mc_pm_md(a.ospa_loc_mean); ...
    'OSPA card (m)',     @(a) mc_pm_md(a.ospa_card_mean); ...
    'OSPA assign (m)',   @(a) mc_pm_md(a.ospa_assign_mean); ...
    '基数误差均值',       @(a) mc_pm_md(a.card_err_mean); ...
    '虚警确认比例',       @(a) mc_pm_md(a.false_confirm_ratio,'pct'); ...
    '中位航迹寿命 (帧)',  @(a) mc_pm_md(a.median_track_life,'f1'); ...
    '碎裂次数',           @(a) mc_pm_md(a.fragmentation_count,'f1'); ...
    'ID Switch',         @(a) mc_pm_md(a.id_switches,'f1') };
for r = 1:size(rows,1)
    fprintf(fid, '| %s |', rows{r,1});
    for i = 1:4
        fprintf(fid, ' %s |', rows{r,2}(packs.(keys{i}).agg));
    end
    fprintf(fid, '\n');
end
fclose(fid);
copyfile(fullfile(out_dir,'exp2_summary.md'), fullfile(tex_dir,'exp2_summary.md'));
end

function plot_exp2_figures(out_dir, tex_dir, packs, keys, names, t_sec)
cols = [0.55 0.55 0.55; 0.85 0.45 0.15; 0.60 0.20 0.55; 0.15 0.35 0.70];

fig1 = figure('Color','w','Position',[60 60 920 400]);
for i=1:4
    mc_plot_series_ci(t_sec, packs.(keys{i}).agg.ospa_k, cols(i,:), 25); hold on;
end
grid on; box on;
xlabel('Time (s)'); ylabel('OSPA (m)');
legend(names, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig1, out_dir, tex_dir, 'fig_exp2_ospa');

fig2 = figure('Color','w','Position',[60 60 920 400]);
plot(t_sec, packs.(keys{1}).agg.n_true_k, 'k-', 'LineWidth', 1.6); hold on;
for i=1:4
    mc_plot_series_ci(t_sec, packs.(keys{i}).agg.n_est_k, cols(i,:), 25);
end
grid on; box on;
xlabel('Time (s)'); ylabel('Number of tracks');
legend(['N_{true}', names], 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig2, out_dir, tex_dir, 'fig_exp2_cardinality');

fig3 = figure('Color','w','Position',[60 60 560 400]);
slist = cell(1,4);
for i=1:4, slist{i} = scale_pct(packs.(keys{i}).agg.false_confirm_ratio); end
mc_bar_with_ci(slist, cols, names);
ylabel('False confirmed track ratio (\%)');
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
export_both(fig3, out_dir, tex_dir, 'fig_exp2_false_confirm');

fig4 = figure('Color','w','Position',[60 60 920 400]);
subplot(1,2,1);
ll = cell(1,4); for i=1:4, ll{i} = packs.(keys{i}).agg.median_track_life; end
mc_bar_with_ci(ll, cols, names);
ylabel('Median track life (frames)'); grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',10);
subplot(1,2,2);
ff = cell(1,4); for i=1:4, ff{i} = packs.(keys{i}).agg.fragmentation_count; end
mc_bar_with_ci(ff, cols, names);
ylabel('Fragmentation count'); grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',10);
export_both(fig4, out_dir, tex_dir, 'fig_exp2_life_frag');
end

function write_exp2_tex_table(tex_dir, packs, keys, n_mc, is_debug)
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
A = packs.(keys{1}).agg; P = packs.(keys{2}).agg;
M = packs.(keys{3}).agg; F = packs.(keys{4}).agg;
fid = fopen(fullfile(tex_dir, 'tab_exp2_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp2 MC table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{实验二：声学观测形式对比（$\\gamma=0.3$，%s）}\n', cap_note);
fprintf(fid, '  \\label{tab:exp2_metrics}\n');
fprintf(fid, '  \\begin{tabular}{lcccc}\n    \\toprule\n');
fprintf(fid, '    指标 & Baseline-R & Point-A & Point-A-Multi & Field-A \\\\\n    \\midrule\n');
fprintf(fid, '    OSPA 均值 (m) & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_mean), mc_pm_tex(P.ospa_mean), mc_pm_tex(M.ospa_mean), mc_pm_tex(F.ospa_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{loc}}$ (m) & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_loc_mean), mc_pm_tex(P.ospa_loc_mean), mc_pm_tex(M.ospa_loc_mean), mc_pm_tex(F.ospa_loc_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{card}}$ (m) & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_card_mean), mc_pm_tex(P.ospa_card_mean), mc_pm_tex(M.ospa_card_mean), mc_pm_tex(F.ospa_card_mean));
fprintf(fid, '    $\\mathrm{OSPA}_{\\mathrm{assign}}$ (m) & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.ospa_assign_mean), mc_pm_tex(P.ospa_assign_mean), mc_pm_tex(M.ospa_assign_mean), mc_pm_tex(F.ospa_assign_mean));
fprintf(fid, '    基数误差均值 & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.card_err_mean), mc_pm_tex(P.card_err_mean), mc_pm_tex(M.card_err_mean), mc_pm_tex(F.card_err_mean));
fprintf(fid, '    虚警确认比例 & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.false_confirm_ratio,'pct'), mc_pm_tex(P.false_confirm_ratio,'pct'), ...
    mc_pm_tex(M.false_confirm_ratio,'pct'), mc_pm_tex(F.false_confirm_ratio,'pct'));
fprintf(fid, '    中位航迹寿命 (帧) & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.median_track_life,'f1'), mc_pm_tex(P.median_track_life,'f1'), ...
    mc_pm_tex(M.median_track_life,'f1'), mc_pm_tex(F.median_track_life,'f1'));
fprintf(fid, '    碎裂次数 & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.fragmentation_count,'f1'), mc_pm_tex(P.fragmentation_count,'f1'), ...
    mc_pm_tex(M.fragmentation_count,'f1'), mc_pm_tex(F.fragmentation_count,'f1'));
fprintf(fid, '    ID Switch & %s & %s & %s & %s \\\\\n', ...
    mc_pm_tex(A.id_switches,'f1'), mc_pm_tex(P.id_switches,'f1'), ...
    mc_pm_tex(M.id_switches,'f1'), mc_pm_tex(F.id_switches,'f1'));
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
