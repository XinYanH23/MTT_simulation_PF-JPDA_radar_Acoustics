%% =========================================================================
%  run_exp5_sparse_acoustic.m
%  实验五：声学站稀疏部署（25 组独立种子；rng(42) 仅调试）
%  站点子集用固定布局种子 LAYOUT_SEED，不随跟踪种子变化。
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

out_dir = fullfile(this_dir, 'experiments', 'exp5_sparse_acoustic');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp5');
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
LAYOUT_SEED = 2024;   % 站点子集固定，与跟踪随机性分离

fprintf('===== 实验五：加载场景与探测 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time); %#ok<NASGU>
pos_all = data1.acousticPos;

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
src_h = median(hs);

n_list = 20:2:60;
n_list = unique([n_list(:)', 33]);   % 偶数加密 + 原场景 33 站
env_nl = getenv('TRACKER_EXP5_N_LIST');
if ~isempty(env_nl)
    n_list = str2num(env_nl); %#ok<ST2NM>
    n_list = unique(n_list(:)');
end
n_max = max(n_list);
buildings = [];
if isfield(data1, 'buildings'); buildings = data1.buildings; end
[pos_pool, extra_pos] = expand_acoustic_catalog(pos_all, buildings, n_max, LAYOUT_SEED);
fprintf('  原站=%d  补点=%d  扫描 n=%s\n', size(pos_all,1), size(extra_pos,1), mat2str(n_list));

packs = struct();
station_idx = struct();
cache_mat = fullfile(out_dir, 'exp5_results.mat');
if exist(cache_mat, 'file')
    old = load(cache_mat, 'packs');
    if isfield(old, 'packs'); packs = old.packs; end
    fprintf('  已加载缓存 %s\n', cache_mat);
end

fprintf('\n########## Radar-only (n=0)  (%d seeds) ##########\n', n_mc);
if isfield(packs, 'n0') && isfield(packs.n0, 'agg')
    fprintf('  跳过（缓存已有）\n');
else
    opts0 = struct('acoustic_mode', 'none', 'verbose', false, 'method', 'softmax');
    packs.n0 = mc_run_case(allDetections, truth, time, cfg0, [], opts0, seeds);
    packs.n0.n_stations = 0;
    packs.n0.label = 'Radar-only';
end

for ni = 1:numel(n_list)
    n_keep = n_list(ni);
    key = sprintf('n%d', n_keep);
    [pos_n, idx] = select_exp5_stations(pos_all, extra_pos, n_keep, LAYOUT_SEED);
    station_idx.(key) = idx;
    if isfield(packs, key) && isfield(packs.(key), 'agg') ...
            && isfield(packs.(key), 'n_stations') && packs.(key).n_stations == n_keep
        fprintf('\n########## Field n=%d  跳过（缓存已有） ##########\n', n_keep);
        continue
    end
    ac_cfg = acoustic_config_from_scenario(pos_n, src_h);
    fprintf('\n########## Field n=%d  (%d seeds) ##########\n', n_keep, n_mc);
    opts = struct('acoustic_mode', 'field', 'verbose', false, 'method', 'softmax');
    packs.(key) = mc_run_case(allDetections, truth, time, cfg0, ac_cfg, opts, seeds);
    packs.(key).n_stations = n_keep;
    packs.(key).label = sprintf('Field-%d', n_keep);
    save(cache_mat, 'packs', 'station_idx', 'n_list', 'truth', 'time', 'cfg0', ...
        'seeds', 'n_mc', 'is_debug', 'LAYOUT_SEED', 'pos_pool', '-v7.3');
end

save(fullfile(out_dir, 'exp5_results.mat'), ...
    'packs', 'station_idx', 'n_list', 'truth', 'time', 'cfg0', ...
    'seeds', 'n_mc', 'is_debug', 'LAYOUT_SEED', 'pos_pool', '-v7.3');

keys = [{'n0'}, arrayfun(@(n) sprintf('n%d', n), n_list, 'UniformOutput', false)];
csv_p = cell(1, numel(keys));
csv_n = cell(1, numel(keys));
for i = 1:numel(keys)
    csv_p{i} = packs.(keys{i});
    csv_n{i} = keys{i};
end
mc_write_seed_csv(fullfile(out_dir, 'exp5_per_seed.csv'), csv_p, csv_n);

fid = fopen(fullfile(out_dir, 'exp5_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验五：声学站稀疏部署（蒙特卡洛）\n\n');
fprintf(fid, '跟踪种子：%d 组 %s；布局种子 LAYOUT_SEED=%d；rng(42) 仅为调试。\n\n', ...
    n_mc, mat2str(seeds), LAYOUT_SEED);
fprintf(fid, '| 站点数 | OSPA | 虚警确认 | 基数误差 | 中位寿命 | 碎裂 | IDS |\n');
fprintf(fid, '|-------:|------|----------|----------|----------|------|-----|\n');
for i = 1:numel(keys)
    a = packs.(keys{i}).agg;
    ns = packs.(keys{i}).n_stations;
    fprintf(fid, '| %d | %s | %s | %s | %s | %s | %s |\n', ns, ...
        mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
        mc_pm_md(a.card_err_mean), mc_pm_md(a.median_track_life,'f1'), ...
        mc_pm_md(a.fragmentation_count,'f1'), mc_pm_md(a.id_switches,'f1'));
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp5_summary.md'), fullfile(tex_dir, 'exp5_summary.md'));

ns = double([0, n_list(:)']);
s_ospa = cell(1, numel(ns));
s_fa   = cell(1, numel(ns));
s_card = cell(1, numel(ns));
for i = 1:numel(ns)
    a = packs.(sprintf('n%d', ns(i))).agg;
    s_ospa{i} = a.ospa_mean;
    s_fa{i}   = scale_pct(a.false_confirm_ratio);
    s_card{i} = a.card_err_mean;
end
[ns_s, ord] = sort(ns);
s_ospa = s_ospa(ord); s_fa = s_fa(ord); s_card = s_card(ord);

fig1 = figure('Color','w','Position',[80 80 720 400]);
mc_errorbar_xy(ns_s, s_ospa, [0.15 0.35 0.70], 'o');
grid on; box on;
xlabel('Number of acoustic stations'); ylabel('OSPA mean (m)');
xticks_n = unique([0, 20:10:60, 33]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', xticks_n);
export_both(fig1, out_dir, tex_dir, 'fig_exp5_ospa_vs_n');

fig2 = figure('Color','w','Position',[80 80 720 400]);
mc_errorbar_xy(ns_s, s_fa, [0.75 0.25 0.25], 's');
grid on; box on;
xlabel('Number of acoustic stations'); ylabel('False confirmed track ratio (%)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', xticks_n);
export_both(fig2, out_dir, tex_dir, 'fig_exp5_fa_vs_n');

fig3 = figure('Color','w','Position',[80 80 720 400]);
mc_errorbar_xy(ns_s, s_card, [0.2 0.5 0.3], 'd');
grid on; box on;
xlabel('Number of acoustic stations'); ylabel('Mean cardinality error');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'XTick', xticks_n);
export_both(fig3, out_dir, tex_dir, 'fig_exp5_card_vs_n');

n_show = intersect([20 33 46 60], ns_s);
if isempty(n_show); n_show = ns_s(ns_s > 0); n_show = n_show(round(linspace(1, numel(n_show), min(4, numel(n_show))))); end
fig4 = figure('Color','w','Position',[80 80 240*numel(n_show) 320]);
for i = 1:numel(n_show)
    subplot(1, numel(n_show), i);
    n_keep = n_show(i);
    [pos_n, ~] = select_exp5_stations(pos_all, extra_pos, n_keep, LAYOUT_SEED);
    scatter(pos_pool(:,1), pos_pool(:,2), 14, [0.82 0.82 0.82], 'filled'); hold on;
    scatter(pos_n(:,1), pos_n(:,2), 36, [0.15 0.35 0.70], 'filled');
    axis equal tight; grid on; box on;
    title(sprintf('n=%d', n_keep));
    xlabel('x (m)'); ylabel('y (m)');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
end
export_both(fig4, out_dir, tex_dir, 'fig_exp5_layout');

write_exp5_tex_table(tex_dir, packs, ns, n_mc, is_debug, cfg0.pf_update_mode);
fprintf('\n===== 实验五完成（n_mc=%d） =====\n', n_mc);

function write_exp5_tex_table(tex_dir, packs, ns, n_mc, is_debug, pf_mode)
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
fid = fopen(fullfile(tex_dir, 'tab_exp5_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp5 MC table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
if nargin < 6 || isempty(pf_mode); pf_mode = 'legacy_likelihood'; end
mode_note = strrep(pf_mode, '_', '\\_');
fprintf(fid, '  \\caption{实验五：站点数扫描（%s，%s）}\n', mode_note, cap_note);
fprintf(fid, '  \\label{tab:exp5_metrics}\n');
fprintf(fid, '  \\begin{tabular}{ccccccc}\n    \\toprule\n');
fprintf(fid, '    站点数 $n$ & OSPA & 虚警确认 & 基数误差 & 中位寿命 & 碎裂 & ID Switch \\\\\n    \\midrule\n');
ns_s = sort(ns);
for i = 1:numel(ns_s)
    a = packs.(sprintf('n%d', ns_s(i))).agg;
    if ns_s(i) == 0
        lab = '0（Radar-only）';
    else
        lab = sprintf('%d', ns_s(i));
    end
    fprintf(fid, '    %s & %s & %s & %s & %s & %s & %s \\\\\n', lab, ...
        mc_pm_tex(a.ospa_mean), mc_pm_tex(a.false_confirm_ratio,'pct'), ...
        mc_pm_tex(a.card_err_mean), mc_pm_tex(a.median_track_life,'f1'), ...
        mc_pm_tex(a.fragmentation_count,'f1'), mc_pm_tex(a.id_switches,'f1'));
end
fprintf(fid, '    \\bottomrule\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);
end

function [pos_n, idx] = select_exp5_stations(pos_base, extra_pos, n_keep, seed)
M = size(pos_base, 1);
if n_keep <= M
    [pos_n, idx] = subsample_acoustic_stations(pos_base, n_keep, seed);
else
    n_ex = n_keep - M;
    if n_ex > size(extra_pos, 1)
        error('select_exp5_stations:Short', '补点不足：需要 %d，只有 %d。', n_ex, size(extra_pos,1));
    end
    pos_n = [pos_base; extra_pos(1:n_ex, :)];
    idx = (1:n_keep)';
end
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
