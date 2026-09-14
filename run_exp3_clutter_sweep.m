%% =========================================================================
%  run_exp3_clutter_sweep.m
%  实验三：雷达虚警/杂波密度扫描（25 组独立种子；rng(42) 仅调试）
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

out_dir = fullfile(this_dir, 'experiments', 'exp3_clutter_sweep');
tex_dir = fullfile(this_dir, 'report', 'figures', 'ch3_exp3');
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

fprintf('===== 实验三：加载场景与探测 =====\n');
load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;
N = length(time);

base_fa = count_radar_fa(allDetections, truth, time, 15);
fprintf('  基线雷达 FA/frame ≈ %.3f（15 m 门）\n', base_fa.fa_per_frame);

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

mu_list = [0, 1, 3, 5, 10];
modes = { ...
    struct('key','RadarOnly', 'acoustic_mode','none',  'label','Radar-only'), ...
    struct('key','Field',     'acoustic_mode','field', 'label','Radar+Acoustic') };

fprintf('  μ 扫描 = %s   n_mc=%d\n', mat2str(mu_list), n_mc);

metrics_all = struct();
csv_packs = {};
csv_names = {};

for ui = 1:numel(mu_list)
    mu = mu_list(ui);
    fprintf('\n========== μ_extra = %.1f  (%d seeds) ==========\n', mu, n_mc);

    for mi = 1:numel(modes)
        md = modes{mi};
        m_list = cell(1, n_mc);
        fa_vec = zeros(n_mc, 1);
        for si = 1:n_mc
            inj_opts = struct( ...
                'mu_fa_extra', mu, ...
                'map_xlim', [0, 1000], ...
                'map_ylim', [0, 1000], ...
                'R_radar', cfg0.R_radar, ...
                'seed', seeds(si)*10 + ui);
            [dets_mu, ~] = inject_radar_clutter(allDetections, inj_opts);
            fa_mu = count_radar_fa(dets_mu, truth, time, 15);
            fa_vec(si) = fa_mu.fa_per_frame;

            fprintf('    %s  [MC %d/%d] rng(%d)\n', md.label, si, n_mc, seeds(si));
            rng(seeds(si), 'twister');
            opts = struct('acoustic_mode', md.acoustic_mode, ...
                          'verbose', false, 'method', 'softmax');
            out = run_tracker_once(dets_mu, truth, time, cfg0, ac_cfg, opts);
            m_list{si} = eval_exp1_metrics(out.tracks, truth, time, out.cfg, ...
                                           out.frame_log, true);
            m_list{si}.seed = seeds(si);
        end
        tag = sprintf('mu%d_%s', round(mu), md.key);
        pack = struct('seeds', seeds, 'per_seed', {m_list});
        pack.agg = mc_aggregate_metrics(m_list);
        pack.agg.fa_per_frame = mc_scalar_stats(fa_vec);
        pack.agg.mu_extra = mu;
        pack.agg.label = md.label;
        metrics_all.(tag) = pack;
        csv_packs{end+1} = pack; %#ok<AGROW>
        csv_names{end+1} = sprintf('mu%g-%s', mu, md.key); %#ok<AGROW>
        a = pack.agg;
        fprintf('    %s  OSPA %s  FA_conf %s  card %s\n', md.label, ...
            mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
            mc_pm_md(a.card_err_mean));
    end
end

save(fullfile(out_dir, 'exp3_results.mat'), ...
    'metrics_all', 'mu_list', 'base_fa', 'truth', 'time', 'cfg0', ...
    'seeds', 'n_mc', 'is_debug', '-v7.3');
mc_write_seed_csv(fullfile(out_dir, 'exp3_per_seed.csv'), csv_packs, csv_names);

write_exp3_summary(out_dir, tex_dir, metrics_all, mu_list, modes, cfg0, base_fa, seeds, n_mc);
plot_exp3_figures(out_dir, tex_dir, metrics_all, mu_list);
write_exp3_tex_table(tex_dir, metrics_all, mu_list, n_mc, is_debug);

fprintf('\n===== 实验三完成（n_mc=%d） =====\n', n_mc);

function fa = count_radar_fa(allDets, truth, time, gate_m)
N = numel(time); Nt = numel(truth);
n_fa = 0; n_radar = 0;
for k = 1:N
    dets = allDets{k};
    if isempty(dets), continue; end
    if ~iscell(dets), dets = num2cell(dets); end
    gt = zeros(0,2);
    for j = 1:Nt
        if k <= size(truth(j).pos3D,1) && ~any(isnan(truth(j).pos3D(k,1:2)))
            gt(end+1,:) = truth(j).pos3D(k,1:2); %#ok<AGROW>
        end
    end
    for d = 1:numel(dets)
        od = dets{d};
        if isempty(od), continue; end
        stype = '';
        try
            stype = char(od.ObjectAttributes.SensorType);
        catch
            continue
        end
        if ~strcmpi(stype,'Radar'), continue; end
        n_radar = n_radar + 1;
        m = od.Measurement(:);
        xy = m(1:2)';
        if isempty(gt)
            n_fa = n_fa + 1;
        else
            dd = sqrt(sum((gt - xy).^2, 2));
            if min(dd) >= gate_m
                n_fa = n_fa + 1;
            end
        end
    end
end
fa = struct('n_radar', n_radar, 'n_fa', n_fa, ...
    'fa_per_frame', n_fa / max(N,1), 'fa_ratio', n_fa / max(n_radar,1));
end

function write_exp3_summary(out_dir, tex_dir, M, mu_list, modes, cfg, base_fa, seeds, n_mc)
fid = fopen(fullfile(out_dir, 'exp3_summary.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '# 实验三：雷达虚警/杂波密度扫描（蒙特卡洛）\n\n');
fprintf(fid, '| 项目 | 取值 |\n|------|------|\n');
fprintf(fid, '| 基线 FA/frame（15 m） | %.3f |\n', base_fa.fa_per_frame);
fprintf(fid, '| 额外杂波 μ | %s |\n', mat2str(mu_list));
fprintf(fid, '| λ_c / γ | %.1e / %.1f |\n', cfg.lambda_c, cfg.pf_ac_gamma_fixed);
fprintf(fid, '| 种子 | 主实验 %d 组 %s；rng(42) 仅为调试 |\n\n', n_mc, mat2str(seeds));
fprintf(fid, '| μ | 方案 | OSPA | 虚警确认 | 基数误差 | 碎裂 | IDS | 中位寿命 |\n');
fprintf(fid, '|--:|------|------|----------|----------|------|-----|----------|\n');
for ui = 1:numel(mu_list)
    mu = mu_list(ui);
    for mi = 1:numel(modes)
        md = modes{mi};
        a = M.(sprintf('mu%d_%s', round(mu), md.key)).agg;
        fprintf(fid, '| %g | %s | %s | %s | %s | %s | %s | %s |\n', mu, md.label, ...
            mc_pm_md(a.ospa_mean), mc_pm_md(a.false_confirm_ratio,'pct'), ...
            mc_pm_md(a.card_err_mean), mc_pm_md(a.fragmentation_count,'f1'), ...
            mc_pm_md(a.id_switches,'f1'), mc_pm_md(a.median_track_life,'f1'));
    end
end
fclose(fid);
copyfile(fullfile(out_dir, 'exp3_summary.md'), fullfile(tex_dir, 'exp3_summary.md'));
end

function plot_exp3_figures(out_dir, tex_dir, M, mu_list)
colR = [0.75 0.25 0.25];
colF = [0.15 0.35 0.70];
n = numel(mu_list);
sR_ospa = cell(1,n); sF_ospa = cell(1,n);
sR_fa = cell(1,n); sF_fa = cell(1,n);
sR_card = cell(1,n); sF_card = cell(1,n);
for ui = 1:n
    mu = mu_list(ui);
    a = M.(sprintf('mu%d_RadarOnly', round(mu))).agg;
    b = M.(sprintf('mu%d_Field', round(mu))).agg;
    sR_ospa{ui} = a.ospa_mean; sF_ospa{ui} = b.ospa_mean;
    sR_fa{ui} = scale_pct(a.false_confirm_ratio);
    sF_fa{ui} = scale_pct(b.false_confirm_ratio);
    sR_card{ui} = a.card_err_mean; sF_card{ui} = b.card_err_mean;
end

fig1 = figure('Color','w','Position',[80 80 720 420]);
mc_errorbar_xy(mu_list, sR_ospa, colR, 'o'); hold on;
mc_errorbar_xy(mu_list, sF_ospa, colF, 's');
grid on; box on;
xlabel('Extra clutter mean \mu (FA/frame)');
ylabel('OSPA mean (m)');
legend({'Radar-only','Radar+Acoustic'}, 'Location', 'northwest');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig1, out_dir, tex_dir, 'fig_exp3_ospa_vs_mu');

fig2 = figure('Color','w','Position',[80 80 720 420]);
mc_errorbar_xy(mu_list, sR_fa, colR, 'o'); hold on;
mc_errorbar_xy(mu_list, sF_fa, colF, 's');
grid on; box on;
xlabel('Extra clutter mean \mu (FA/frame)');
ylabel('False confirmed track ratio (%)');
legend({'Radar-only','Radar+Acoustic'}, 'Location', 'northwest');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig2, out_dir, tex_dir, 'fig_exp3_fa_vs_mu');

fig3 = figure('Color','w','Position',[80 80 720 420]);
mc_errorbar_xy(mu_list, sR_card, colR, 'o'); hold on;
mc_errorbar_xy(mu_list, sF_card, colF, 's');
grid on; box on;
xlabel('Extra clutter mean \mu (FA/frame)');
ylabel('Mean cardinality error');
legend({'Radar-only','Radar+Acoustic'}, 'Location', 'northwest');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig3, out_dir, tex_dir, 'fig_exp3_card_vs_mu');

d_ospa = cell(1,n); d_fa = cell(1,n);
for ui = 1:n
    mu = mu_list(ui);
    a = M.(sprintf('mu%d_RadarOnly', round(mu)));
    b = M.(sprintf('mu%d_Field', round(mu)));
    d_ospa{ui} = mc_scalar_stats(a.agg.raw.ospa_mean - b.agg.raw.ospa_mean);
    d_fa{ui}   = mc_scalar_stats(100*(a.agg.raw.false_confirm_ratio - b.agg.raw.false_confirm_ratio));
end
fig4 = figure('Color','w','Position',[80 80 720 420]);
yyaxis left;
mc_errorbar_xy(mu_list, d_ospa, [0.2 0.5 0.3], 'd');
ylabel('\Delta OSPA (Radar-only - Field) (m)');
yyaxis right;
mc_errorbar_xy(mu_list, d_fa, [0.6 0.4 0.1], '^');
ylabel('\Delta False-confirm (pp)');
grid on; box on;
xlabel('Extra clutter mean \mu (FA/frame)');
legend({'\Delta OSPA','\Delta FA confirm'}, 'Location', 'northwest');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
export_both(fig4, out_dir, tex_dir, 'fig_exp3_delta_benefit');
end

function write_exp3_tex_table(tex_dir, M, mu_list, n_mc, is_debug)
if is_debug
    cap_note = '调试种子 \texttt{rng(42)}（非主实验）';
else
    cap_note = sprintf('%d 组独立随机种子，单元格为均值 $\\pm$ 标准差', n_mc);
end
fid = fopen(fullfile(tex_dir, 'tab_exp3_metrics.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Auto-generated Exp3 MC table\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{实验三：杂波扫描主要指标（$\\gamma=0.3$，%s）}\n', cap_note);
fprintf(fid, '  \\label{tab:exp3_metrics}\n');
fprintf(fid, '  \\setlength{\\tabcolsep}{3.5pt}\n');
fprintf(fid, '  \\begin{tabular}{clcccc}\n    \\toprule\n');
fprintf(fid, '    $\\mu$ & 方案 & OSPA & 虚警确认 & 基数误差 & 碎裂 \\\\\n    \\midrule\n');
for ui = 1:numel(mu_list)
    mu = mu_list(ui);
    a = M.(sprintf('mu%d_RadarOnly', round(mu))).agg;
    b = M.(sprintf('mu%d_Field', round(mu))).agg;
    fprintf(fid, '    %g & Radar-only & %s & %s & %s & %s \\\\\n', mu, ...
        mc_pm_tex(a.ospa_mean), mc_pm_tex(a.false_confirm_ratio,'pct'), ...
        mc_pm_tex(a.card_err_mean), mc_pm_tex(a.fragmentation_count,'f1'));
    fprintf(fid, '    %g & +Field & %s & %s & %s & %s \\\\\n', mu, ...
        mc_pm_tex(b.ospa_mean), mc_pm_tex(b.false_confirm_ratio,'pct'), ...
        mc_pm_tex(b.card_err_mean), mc_pm_tex(b.fragmentation_count,'f1'));
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
