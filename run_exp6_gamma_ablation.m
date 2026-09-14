%% =========================================================================
%  run_exp6_gamma_ablation.m
%  声学融合权重 γ 消融：固定其余参数，仅扫描 pf_ac_gamma_fixed。
%  主实验 25 组独立种子；rng(42) 仅为调试。
%
%  融合公式（与 PF/KF 共用）：
%    ℓ_a = γ L_k(x) + (1-γ)
%  γ=0 时更新步不吃场（ℓ_a≡1），但场生成与 η_a 起始门限保持不变，
%  从而把“更新权重”从“有无声学场”中隔离出来。
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

out_dir = fullfile(this_dir, 'experiments', 'exp6_gamma_ablation');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp6');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

MC_DEBUG = false;
n_mc     = 25;
if strcmp(getenv('TRACKER_MC_DEBUG'), '1'); MC_DEBUG = true; end
env_n = str2double(getenv('TRACKER_N_MC'));
if ~isnan(env_n) && env_n >= 1; n_mc = env_n; end
[seeds, n_mc, is_debug] = mc_resolve_seeds(MC_DEBUG, n_mc); %#ok<ASGLU>

gamma_list = [0, 0.1, 0.3, 0.5, 0.7, 1.0];
env_g = getenv('TRACKER_GAMMA_LIST');
if ~isempty(env_g)
    gamma_list = str2num(env_g); %#ok<ST2NM>
    gamma_list = gamma_list(:)';
end

fprintf('===== 实验六：加载场景与探测 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time); %#ok<NASGU>

cfg0 = config_tracker();
cfg0.use_pf = true;
cfg0.save_mat = false;
cfg0.jpda_radar_only = true;
cfg0.birth_radar_only = true;
cfg0.pf_update_mode   = 'legacy_likelihood';
cfg0.pf_ac_gamma_mode = 'fixed';

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));
fprintf('  UAV=%d  声学站=%d  γ 扫描=%s\n', numel(truth), ...
    size(data1.acousticPos,1), mat2str(gamma_list));

packs = struct();
keys  = cell(1, numel(gamma_list));
opts  = struct('acoustic_mode', 'field', 'verbose', false, 'method', 'softmax');
for gi = 1:numel(gamma_list)
    g = gamma_list(gi);
    key = gamma_key(g);
    keys{gi} = key;
    fprintf('\n########## γ=%.2f  (%d seeds) ##########\n', g, n_mc);
    cfg = cfg0;
    cfg.pf_ac_gamma_fixed = g;
    packs.(key) = mc_run_case(allDetections, truth, time, cfg, ac_cfg, opts, seeds);
    packs.(key).gamma = g;
end

save(fullfile(out_dir, 'exp6_results.mat'), ...
    'packs', 'gamma_list', 'keys', 'truth', 'time', 'cfg0', ...
    'seeds', 'n_mc', 'is_debug', '-v7.3');

csv_p = cell(1, numel(keys));
csv_n = cell(1, numel(keys));
for i = 1:numel(keys)
    csv_p{i} = packs.(keys{i});
    csv_n{i} = sprintf('gamma_%.2f', packs.(keys{i}).gamma);
end
mc_write_seed_csv(fullfile(out_dir, 'exp6_per_seed.csv'), csv_p, csv_n);

fid = fopen(fullfile(out_dir, 'exp6_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验六：声学融合权重 γ 消融（蒙特卡洛）\n\n');
fprintf(fid, '种子：主实验 %d 组 %s；rng(42) 仅为调试。\n', n_mc, mat2str(seeds));
fprintf(fid, '固定：Field-A / PF / softmax / η_a / 航迹管理；仅扫描 γ。\n\n');
fprintf(fid, '| γ | OSPA | 虚警确认 | 平均寿命 | 基数误差 |\n');
fprintf(fid, '|--:|------|----------|----------|----------|\n');
for i = 1:numel(keys)
    a = packs.(keys{i}).agg;
    fprintf(fid, '| %.2f | %s | %s | %s | %s |\n', packs.(keys{i}).gamma, ...
        mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
        mc_pm_md(a.mean_track_life,'f1'), mc_pm_md(a.card_err_mean));
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp6_summary.md'), fullfile(tex_dir, 'exp6_summary.md'));

g_x = gamma_list;
s_ospa = cell(1, numel(keys));
s_fa   = cell(1, numel(keys));
s_life = cell(1, numel(keys));
s_card = cell(1, numel(keys));
for i = 1:numel(keys)
    a = packs.(keys{i}).agg;
    s_ospa{i} = a.ospa_mean;
    s_fa{i}   = scale_pct(a.false_confirm_ratio);
    s_life{i} = a.mean_track_life;
    s_card{i} = a.card_err_mean;
end

fig = figure('Color','w','Position',[60 60 980 720]);
subplot(2,2,1);
mc_errorbar_xy(g_x, s_ospa, [0.15 0.35 0.70], 'o'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('OSPA mean (m)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);

subplot(2,2,2);
mc_errorbar_xy(g_x, s_fa, [0.75 0.25 0.25], 's'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('False confirmed track ratio (%)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);

subplot(2,2,3);
mc_errorbar_xy(g_x, s_life, [0.20 0.50 0.30], 'd'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('Mean track life (frames)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);

subplot(2,2,4);
mc_errorbar_xy(g_x, s_card, [0.55 0.35 0.15], '^'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('Mean cardinality error');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);
export_both(fig, out_dir, tex_dir, 'fig_exp6_gamma_sweep');

fig1 = figure('Color','w','Position',[80 80 720 380]);
mc_errorbar_xy(g_x, s_ospa, [0.15 0.35 0.70], 'o'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('OSPA mean (m)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);
export_both(fig1, out_dir, tex_dir, 'fig_exp6_ospa');

fig2 = figure('Color','w','Position',[80 80 720 380]);
mc_errorbar_xy(g_x, s_fa, [0.75 0.25 0.25], 's'); hold on;
xline(0.3, 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('\gamma'); ylabel('False confirmed track ratio (%)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', g_x);
export_both(fig2, out_dir, tex_dir, 'fig_exp6_fa');

write_exp6_tex_table(tex_dir, packs, keys, n_mc, is_debug);
fprintf('\n===== 实验六完成（n_mc=%d, n_gamma=%d） =====\n', n_mc, numel(gamma_list));

function key = gamma_key(g)
key = sprintf('g%d', round(g * 100));
end

function write_exp6_tex_table(tex_dir, packs, keys, n_mc, is_debug)
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
fid = fopen(fullfile(tex_dir, 'tab_exp6_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp6 gamma ablation table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, ['  \\caption{声学融合权重 $\\gamma$ 消融（Field-A / PF，其余参数固定；%s。' ...
    '虚线对应本文选用的 $\\gamma=0.3$）}\n'], cap_note);
fprintf(fid, '  \\label{tab:exp6_gamma}\n');
fprintf(fid, '  \\begin{tabular}{ccccc}\n    \\toprule\n');
fprintf(fid, '    $\\gamma$ & OSPA (m) & 虚警确认 & 平均航迹寿命 (帧) & 基数误差 \\\\\n    \\midrule\n');
for i = 1:numel(keys)
    g = packs.(keys{i}).gamma;
    a = packs.(keys{i}).agg;
    glab = sprintf('%.1f', g);
    if abs(g - 0.3) < 1e-9
        glab = [char(92) 'textbf{0.3}'];
    end
    fprintf(fid, '    %s & %s & %s & %s & %s \\\\\n', glab, ...
        mc_pm_tex(a.ospa_mean), mc_pm_tex(a.false_confirm_ratio,'pct'), ...
        mc_pm_tex(a.mean_track_life,'f1'), mc_pm_tex(a.card_err_mean));
end
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
