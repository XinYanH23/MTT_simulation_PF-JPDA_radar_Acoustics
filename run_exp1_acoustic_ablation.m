%% =========================================================================
%  run_exp1_acoustic_ablation.m
%  实验一：有无声学概率观测对比（Radar-only vs Radar+Acoustic Field）
%
%  主实验：25 组独立互不重复随机种子（mc_experiment_seeds），其余参数固定。
%  调试：将 MC_DEBUG=true，仅 rng(42) 跑一次；rng(42) 不是主实验结果。
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

out_dir = fullfile(this_dir, 'experiments', 'exp1_acoustic_ablation');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp1');
env_out = getenv('TRACKER_OUT_DIR');
env_tex = getenv('TRACKER_TEX_DIR');
if ~isempty(env_out); out_dir = env_out; end
if ~isempty(env_tex); tex_dir = env_tex; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

%% ── 蒙特卡洛种子 ────────────────────────────────────────────────────────
MC_DEBUG = false;     % true → rng(42) 单次调试；主实验必须为 false
n_mc     = 25;        % 主实验组数（20–30，论文取 25）
if strcmp(getenv('TRACKER_MC_DEBUG'), '1'); MC_DEBUG = true; end
env_n = str2double(getenv('TRACKER_N_MC'));
if ~isnan(env_n) && env_n >= 1; n_mc = env_n; end
[seeds, n_mc, is_debug] = mc_resolve_seeds(MC_DEBUG, n_mc); %#ok<ASGLU>

%% ── 加载数据 ─────────────────────────────────────────────────────────────
fprintf('===== 实验一：加载场景与探测 =====\n');
s2 = fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat');
s1 = fullfile(this_dir, 'Scenario', 'Step1_Data.mat');
assert(exist(s2,'file')==2 && exist(s1,'file')==2, '缺少 Step1/Step2 数据');
load(s2, 'allDetections', 'time');
data1 = load(s1);
truth = data1.truth;
N = length(time);
fprintf('  帧数=%d  UAV=%d\n', N, numel(truth));

%% ── 公共配置（两方案完全一致，不随种子改变）──────────────────────────
cfg = config_tracker();
cfg.use_pf = true;
cfg.save_mat = false;
cfg.jpda_radar_only = true;
cfg.birth_radar_only = true;
cfg.pf_update_mode = 'legacy_likelihood';
env_mode = getenv('TRACKER_PF_UPDATE_MODE');
if ~isempty(env_mode); cfg.pf_update_mode = env_mode; end
fprintf('  pf_update_mode=%s\n  out_dir=%s\n', cfg.pf_update_mode, out_dir);

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
src_h = median(hs);
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, src_h);
fprintf('  声学站点数=%d  src_h=%.1f m\n', size(data1.acousticPos,1), src_h);

%% ── 方案 A：Radar-only ──────────────────────────────────────────────────
fprintf('\n########## 方案 A：Radar-only  (%d seeds) ##########\n', n_mc);
cfg_A = cfg;
cfg_A.birth_acoustic_eta = 0;
opts_A = struct('enable_acoustic', false, 'verbose', false, 'method', 'softmax');
mc_A = mc_run_case(allDetections, truth, time, cfg_A, ac_cfg, opts_A, seeds);

%% ── 方案 B：Radar + Acoustic Field ──────────────────────────────────────
fprintf('\n########## 方案 B：Radar+Acoustic  (%d seeds) ##########\n', n_mc);
cfg_B = cfg;
opts_B = struct('enable_acoustic', true, 'verbose', false, 'method', 'softmax');
mc_B = mc_run_case(allDetections, truth, time, cfg_B, ac_cfg, opts_B, seeds);

%% ── 保存 ────────────────────────────────────────────────────────────────
t_sec = time(:)';
if max(t_sec) < 1; t_sec = (0:N-1) * cfg.dt; end

summary = struct();
summary.A = mc_A.agg;
summary.B = mc_B.agg;
summary.seeds = seeds;
summary.n_mc = n_mc;
summary.is_debug = is_debug;
summary.cfg = cfg;
summary.ac_cfg_info = struct( ...
    'n_stations', size(data1.acousticPos,1), ...
    'A', ac_cfg.A(1), 'n', ac_cfg.n(1), 'nf', ac_cfg.noise_floor_db(1), ...
    'src_h', src_h, 'det_threshold', ac_cfg.det_threshold);

save(fullfile(out_dir, 'exp1_results.mat'), ...
    'mc_A', 'mc_B', 'summary', 'truth', 'time', 'cfg', 'seeds', 'n_mc', 'is_debug', '-v7.3');
mc_write_seed_csv(fullfile(out_dir, 'exp1_per_seed.csv'), {mc_A, mc_B}, ...
    {'Radar-only', 'Radar+Acoustic'});
fprintf('结果已保存: %s\n', fullfile(out_dir, 'exp1_results.mat'));

print_mc_headline('Radar-only', mc_A.agg);
print_mc_headline('Radar+Acoustic', mc_B.agg);

%% ── 图 / 表 / 报告 ──────────────────────────────────────────────────────
plot_exp1_mc(out_dir, tex_dir, mc_A.agg, mc_B.agg, t_sec);
write_exp1_report(out_dir, tex_dir, summary, cfg, ac_cfg, data1);
write_exp1_tex_table(tex_dir, mc_A.agg, mc_B.agg, n_mc, is_debug, cfg.pf_update_mode);

fprintf('\n===== 实验一完成（n_mc=%d） =====\n', n_mc);
fprintf('图表与报告目录: %s\n', out_dir);

%% ===== helpers =====
function print_mc_headline(name, a)
fprintf('  %s  OSPA %s  card %s  FA %s  frag %s  IDS %s  life(med) %s\n', ...
    name, mc_pm_md(a.ospa_mean), mc_pm_md(a.card_err_mean), ...
    mc_pm_md(a.false_confirm_ratio, 'pct'), ...
    mc_pm_md(a.fragmentation_count, 'f1'), ...
    mc_pm_md(a.id_switches, 'f1'), ...
    mc_pm_md(a.median_track_life, 'f1'));
end

function plot_exp1_mc(out_dir, tex_dir, A, B, t_sec)
colA = [0.75 0.25 0.25];
colB = [0.15 0.35 0.70];

fig1 = figure('Color','w','Position',[80 80 900 420]);
mc_plot_series_ci(t_sec, A.ospa_k, colA, 25); hold on;
mc_plot_series_ci(t_sec, B.ospa_k, colB, 25);
grid on; box on;
xlabel('Time (s)'); ylabel('OSPA (m)');
title(sprintf('OSPA mean \\pm 95%% CI  (n_{mc}=%d)', A.n_mc));
legend({sprintf('Radar-only (%s)', mc_pm_md(A.ospa_mean)), ...
        sprintf('Radar+Acoustic (%s)', mc_pm_md(B.ospa_mean))}, ...
        'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig1, out_dir, tex_dir, 'fig_exp1_ospa');

fig2 = figure('Color','w','Position',[80 80 900 420]);
plot(t_sec, A.n_true_k, 'k-', 'LineWidth', 1.6); hold on;
mc_plot_series_ci(t_sec, A.n_est_k, colA, 25);
mc_plot_series_ci(t_sec, B.n_est_k, colB, 25);
grid on; box on;
xlabel('Time (s)'); ylabel('Number of tracks');
legend({'$N_{\mathrm{true}}$', 'Radar-only $\hat{N}$', 'Radar+Acoustic $\hat{N}$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig2, out_dir, tex_dir, 'fig_exp1_cardinality');

fig3 = figure('Color','w','Position',[80 80 520 420]);
sA = A.false_confirm_ratio; sA.mean = 100*sA.mean; sA.std = 100*sA.std;
sA.ci_lo = 100*sA.ci_lo; sA.ci_hi = 100*sA.ci_hi;
sB = B.false_confirm_ratio; sB.mean = 100*sB.mean; sB.std = 100*sB.std;
sB.ci_lo = 100*sB.ci_lo; sB.ci_hi = 100*sB.ci_hi;
mc_bar_with_ci({sA, sB}, [colA; colB], {'Radar-only', 'Radar+Acoustic'});
ylabel('False confirmed track ratio (\%)');
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig3, out_dir, tex_dir, 'fig_exp1_false_confirm');

fig4 = figure('Color','w','Position',[80 80 900 420]);
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
export_both(fig4, out_dir, tex_dir, 'fig_exp1_track_life');
end

function write_exp1_report(out_dir, tex_dir, summary, cfg, ac_cfg, data1)
A = summary.A; B = summary.B;
fid = fopen(fullfile(out_dir, 'exp1_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验一：有无声学概率观测对比（蒙特卡洛）\n\n');
fprintf(fid, '## 1. 参数设置\n\n');
fprintf(fid, '| 项目 | 取值 |\n|------|------|\n');
fprintf(fid, '| 随机种子 | 主实验 %d 组独立种子 %s；rng(42) 仅为调试 |\n', ...
    summary.n_mc, mat2str(summary.seeds));
fprintf(fid, '| 场景 | Urban 1000m×1000m，5 UAV |\n');
fprintf(fid, '| 雷达 | 2 站，JPDA 仅雷达点迹 |\n');
fprintf(fid, '| 声学站 | %d（方案B启用） |\n', size(data1.acousticPos,1));
fprintf(fid, '| IMM | 双 DWNA，q=(%.1f, %.1f)，Π 默认 |\n', cfg.cv_q, cfg.ca_q);
fprintf(fid, '| PF | N=%d，ESS 阈值 %.1f |\n', cfg.pf_N, cfg.pf_resample_thresh);
fprintf(fid, '| JPDA | Pd=%.2f，λ_c=%.1e，gate χ²=%.3f |\n', cfg.Pd, cfg.lambda_c, cfg.gate_chi2);
fprintf(fid, '| M/N | %d/%d |\n', cfg.mn_M, cfg.mn_N);
fprintf(fid, '| del_confirmed | %d |\n', cfg.del_confirmed);
fprintf(fid, '| birth_suppress_radius | %.0f m |\n', cfg.birth_suppress_radius);
fprintf(fid, '| 方案A | enable_acoustic=false，η_a=0 |\n');
fprintf(fid, '| pf_update_mode | %s |\n', cfg.pf_update_mode);
fprintf(fid, '| 方案B | enable_acoustic=true，η_a=%.0e，softmax 场 |\n', cfg.birth_acoustic_eta);
fprintf(fid, '| 声学传播 | A=%.0f，n=%.0f，nf=%.0f（与 Step2 对齐） |\n\n', ...
    ac_cfg.A(1), ac_cfg.n(1), ac_cfg.noise_floor_db(1));

fprintf(fid, '## 2. 实验结果（均值 ± 标准差，n=%d）\n\n', summary.n_mc);
fprintf(fid, '| 指标 | Radar-only (A) | Radar+Acoustic (B) |\n');
fprintf(fid, '|------|----------------|--------------------|\n');
fprintf(fid, '| OSPA 均值 (m) | %s | %s |\n', mc_pm_md(A.ospa_mean), mc_pm_md(B.ospa_mean));
fprintf(fid, '| OSPA loc (m) | %s | %s |\n', mc_pm_md(A.ospa_loc_mean), mc_pm_md(B.ospa_loc_mean));
fprintf(fid, '| OSPA card (m) | %s | %s |\n', mc_pm_md(A.ospa_card_mean), mc_pm_md(B.ospa_card_mean));
fprintf(fid, '| OSPA assign (m) | %s | %s |\n', mc_pm_md(A.ospa_assign_mean), mc_pm_md(B.ospa_assign_mean));
fprintf(fid, '| OSPA RMSE (m) | %s | %s |\n', mc_pm_md(A.ospa_rmse), mc_pm_md(B.ospa_rmse));
fprintf(fid, '| 基数误差均值 | %s | %s |\n', mc_pm_md(A.card_err_mean), mc_pm_md(B.card_err_mean));
fprintf(fid, '| 峰值确认航迹数 | %s | %s |\n', mc_pm_md(A.peak_confirmed, 'f1'), mc_pm_md(B.peak_confirmed, 'f1'));
fprintf(fid, '| 虚警确认比例 | %s | %s |\n', mc_pm_md(A.false_confirm_ratio, 'pct'), mc_pm_md(B.false_confirm_ratio, 'pct'));
fprintf(fid, '| 平均航迹寿命 (帧) | %s | %s |\n', mc_pm_md(A.mean_track_life, 'f1'), mc_pm_md(B.mean_track_life, 'f1'));
fprintf(fid, '| 中位航迹寿命 (帧) | %s | %s |\n', mc_pm_md(A.median_track_life, 'f1'), mc_pm_md(B.median_track_life, 'f1'));
fprintf(fid, '| 碎裂次数 | %s | %s |\n', mc_pm_md(A.fragmentation_count, 'f1'), mc_pm_md(B.fragmentation_count, 'f1'));
fprintf(fid, '| ID Switch | %s | %s |\n', mc_pm_md(A.id_switches, 'f1'), mc_pm_md(B.id_switches, 'f1'));
fprintf(fid, '| 位置 RMSE (m，附带) | %s | %s |\n\n', mc_pm_md(A.pos_rmse), mc_pm_md(B.pos_rmse));
fprintf(fid, '逐种子明细见 `exp1_per_seed.csv`。95%% CI 见图中误差棒。\n');
fclose(fid);
copyfile(fullfile(out_dir, 'exp1_summary.md'), fullfile(tex_dir, 'exp1_summary.md'));
end

function write_exp1_tex_table(tex_dir, A, B, n_mc, is_debug, pf_mode)
if nargin < 6 || isempty(pf_mode); pf_mode = 'legacy_likelihood'; end
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
if strcmp(pf_mode, 'jpda_posterior')
    cap_method = 'JPDA 后验 $\\times$ 声学相对似然';
else
    cap_method = '$\\gamma=0.3$';
end
fid = fopen(fullfile(tex_dir, 'tab_exp1_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated by run_exp1_acoustic_ablation.m\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{实验一：Radar-only 与 Radar+Acoustic Field 对比（%s，%s）}\n', cap_method, cap_note);
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
fprintf(fid, '    峰值确认航迹数 & %s & %s \\\\\n', mc_pm_tex(A.peak_confirmed, 'f1'), mc_pm_tex(B.peak_confirmed, 'f1'));
fprintf(fid, '    虚警确认比例 & %s & %s \\\\\n', mc_pm_tex(A.false_confirm_ratio, 'pct'), mc_pm_tex(B.false_confirm_ratio, 'pct'));
fprintf(fid, '    中位航迹寿命 (帧) & %s & %s \\\\\n', mc_pm_tex(A.median_track_life, 'f1'), mc_pm_tex(B.median_track_life, 'f1'));
fprintf(fid, '    碎裂次数 & %s & %s \\\\\n', mc_pm_tex(A.fragmentation_count, 'f1'), mc_pm_tex(B.fragmentation_count, 'f1'));
fprintf(fid, '    ID Switch & %s & %s \\\\\n', mc_pm_tex(A.id_switches, 'f1'), mc_pm_tex(B.id_switches, 'f1'));
fprintf(fid, '    位置 RMSE (m，附带) & %s & %s \\\\\n', mc_pm_tex(A.pos_rmse), mc_pm_tex(B.pos_rmse));
fprintf(fid, '    \\bottomrule\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);
end

function export_both(fig, out_dir, tex_dir, name)
exportgraphics(fig, fullfile(out_dir, [name '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, [name '.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(tex_dir, [name '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(tex_dir, [name '.pdf']), 'ContentType', 'vector');
fprintf('  saved %s\n', name);
end
