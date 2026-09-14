%% 实验一：各 MC 种子的确认航迹数 + 均值/众数
% 这里的“粒子”是 25 次独立蒙特卡洛（种子），不是 PF 内部 500 个粒子。
% 航迹条数由 track_manage 逐帧给出，每个种子每帧都是整数。
clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'Evaluation'));

mat_b = fullfile(this_dir, 'methods_v2_bayes', 'experiments', ...
    'exp1_acoustic_ablation', 'exp1_results.mat');
assert(exist(mat_b,'file')==2, '缺少方法 B 实验一 mat');
S = load(mat_b, 'mc_A', 'mc_B', 'time', 'n_mc');
t = S.time(:)';
if max(t) < 1; t = (0:numel(t)-1) * 0.1; end

out_dir = fullfile(this_dir, 'methods_v2_bayes', 'experiments', 'exp1_acoustic_ablation');
tex_dir = fullfile(this_dir, 'methods_v2_bayes', 'figures', 'ch3_exp1');
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

cases = { ...
    struct('pack', S.mc_A, 'name', 'Radar-only', 'file', 'radar'), ...
    struct('pack', S.mc_B, 'name', 'Radar+Field', 'file', 'field') };

for ci = 1:numel(cases)
    cs = cases{ci};
    [nest, ntrue, mu, md] = nest_from_pack(cs.pack);
    n_mc = size(nest, 1);

    ymax = max([nest(:); ntrue(:)]) + 1.2;
    cmap = lines(n_mc);

    fig = figure('Color','w','Position',[60 60 980 500]);
    hold on;
    for i = 1:n_mc
        plot(t, nest(i,:), 'Color', [cmap(i,:) 0.55], 'LineWidth', 1.0);
    end
    hT = plot(t, ntrue, 'k-', 'LineWidth', 2.4);
    hM = plot(t, mu, 'Color', [0.80 0.15 0.12], 'LineWidth', 2.0);
    hD = plot(t, md, 'Color', [0.05 0.40 0.18], 'LineWidth', 2.2, 'LineStyle', '-.');
    grid on; box on;
    xlabel('Time (s)'); ylabel('Number of confirmed tracks');
    title(sprintf('%s  (method B, %d MC seeds; each thin line = one seed)', cs.name, n_mc));
    legend([hT hM hD], {'N true', 'mean over seeds', 'mode over seeds'}, ...
        'Location', 'northeast', 'Interpreter', 'none');
    ylim([-0.2, ymax]);
    set(gca, 'FontSize', 11, 'YTick', 0:1:20);
    export_png(fig, out_dir, tex_dir, ['fig_exp1_N_seeds_' cs.file]);

    fig2 = figure('Color','w','Position',[60 60 980 440]);
    plot(t, ntrue, 'k-', 'LineWidth', 2.4); hold on;
    plot(t, mu, 'Color', [0.80 0.15 0.12], 'LineWidth', 1.8);
    plot(t, md, 'Color', [0.05 0.40 0.18], 'LineWidth', 2.2);
    grid on; box on;
    xlabel('Time (s)'); ylabel('Number of confirmed tracks');
    title(sprintf('%s: mean vs mode of confirmed-track count', cs.name));
    legend({'N true', 'mean (can be fractional)', 'mode (integer)'}, ...
        'Location', 'northeast', 'Interpreter', 'none');
    ylim([-0.2, max([mu(:); md(:); ntrue(:)]) + 1.2]);
    set(gca, 'FontSize', 11, 'YTick', 0:1:12);
    export_png(fig2, out_dir, tex_dir, ['fig_exp1_N_mode_' cs.file]);

    nmax = max(nest(:));
    heat = zeros(nmax + 1, numel(t));
    for k = 1:numel(t)
        heat(:, k) = histcounts(nest(:, k), -0.5:1:(nmax + 0.5))';
    end
    fig3 = figure('Color','w','Position',[60 60 980 360]);
    imagesc(t, 0:nmax, heat);
    axis xy; colormap(parula); cb = colorbar;
    cb.Label.String = 'number of seeds';
    hold on;
    plot(t, ntrue, 'w-', 'LineWidth', 2.0);
    xlabel('Time (s)'); ylabel('Confirmed-track count N');
    title(sprintf('%s: how many seeds report each integer N(t)', cs.name));
    set(gca, 'FontSize', 11, 'YTick', 0:nmax);
    export_png(fig3, out_dir, tex_dir, ['fig_exp1_N_heat_' cs.file]);

    agree = mean(md == ntrue) * 100;
    mae_mu = mean(abs(mu - ntrue));
    mae_md = mean(abs(double(md) - ntrue));
    fprintf('%s  mode==N_true %.1f%%  MAE mean %.3f  MAE mode %.3f\n', ...
        cs.name, agree, mae_mu, mae_md);
end

% 方法 A 生命周期 25 种子（正式 A 实验一 mat 已坏，用这份看 A 的细线）
mat_a = fullfile(this_dir, 'experiments', 'exp1_ac_lifecycle', 'exp1_ac_lifecycle.mat');
if exist(mat_a, 'file')
    A = load(mat_a, 'packs', 'time');
    if isfield(A, 'packs') && isfield(A.packs, 'A_life')
        tt = A.time(:)';
        if max(tt) < 1; tt = (0:numel(tt)-1) * 0.1; end
        [nest, ntrue, mu, md] = nest_from_pack(A.packs.A_life);
        n_mc = size(nest, 1);
        cmap = lines(n_mc);
        fig = figure('Color','w','Position',[60 60 980 500]);
        hold on;
        for i = 1:n_mc
            plot(tt, nest(i,:), 'Color', [cmap(i,:) 0.55], 'LineWidth', 1.0);
        end
        hT = plot(tt, ntrue, 'k-', 'LineWidth', 2.4);
        hM = plot(tt, mu, 'Color', [0.80 0.15 0.12], 'LineWidth', 2.0);
        hD = plot(tt, md, 'Color', [0.05 0.40 0.18], 'LineWidth', 2.2, 'LineStyle', '-.');
        grid on; box on;
        xlabel('Time (s)'); ylabel('Number of confirmed tracks');
        title(sprintf('Method A + lifecycle  (%d MC seeds; each thin line = one seed)', n_mc));
        legend([hT hM hD], {'N true', 'mean over seeds', 'mode over seeds'}, ...
            'Location', 'northeast', 'Interpreter', 'none');
        set(gca, 'FontSize', 11);
        out_a = fullfile(this_dir, 'experiments', 'exp1_ac_lifecycle');
        export_png(fig, out_a, out_a, 'fig_exp1_N_seeds_A_life');
        fprintf('A+life  mode==N_true %.1f%%  MAE mean %.3f  MAE mode %.3f\n', ...
            mean(md==ntrue)*100, mean(abs(mu-ntrue)), mean(abs(double(md)-ntrue)));
    end
end

function [nest, ntrue, mu, md] = nest_from_pack(pack)
ps = pack.per_seed;
n = numel(ps);
N = numel(ps{1}.n_est_k);
nest = zeros(n, N);
for i = 1:n
    nest(i, :) = ps{i}.n_est_k;
end
ntrue = pack.agg.n_true_k(:)';
if numel(ntrue) ~= N
    ntrue = ps{1}.n_true_k(:)';
end
mu = mean(nest, 1);
md = zeros(1, N);
for k = 1:N
    md(k) = mode(nest(:, k));
end
end

function export_png(fig, out_dir, tex_dir, name)
exportgraphics(fig, fullfile(out_dir, [name '.png']), 'Resolution', 200);
exportgraphics(fig, fullfile(tex_dir, [name '.png']), 'Resolution', 200);
fprintf('saved %s\n', name);
end
