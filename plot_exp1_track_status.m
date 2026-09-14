%% 单次仿真：待定 / 确认 / GT / 删除 / 待定+确认-删除
% 实验一 MC 结果只存了确认数，这里重跑方法 B Field、种子 1008。
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

s2 = fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat');
s1 = fullfile(this_dir, 'Scenario', 'Step1_Data.mat');
load(s2, 'allDetections', 'time');
data1 = load(s1);
truth = data1.truth;
t = time(:)';
if max(t) < 1; t = (0:numel(t)-1) * 0.1; end
N = numel(t);

cfg = config_tracker();
cfg.use_pf = true;
cfg.save_mat = false;
cfg.jpda_radar_only = true;
cfg.birth_radar_only = true;
cfg.pf_update_mode = 'jpda_posterior';

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));

seed = 1008;
rng(seed, 'twister');
opts = struct('enable_acoustic', true, 'verbose', false, 'method', 'softmax');
fprintf('Running method B Field  seed=%d ...\n', seed);
out = run_tracker_once(allDetections, truth, time, cfg, ac_cfg, opts);
fl = out.frame_log;

n_true = zeros(1, N);
for k = 1:N
    for j = 1:numel(truth)
        raw = truth(j).pos3D;
        if k <= size(raw,1) && ~any(isnan(raw(k,1:2)))
            n_true(k) = n_true(k) + 1;
        end
    end
end

n_t = fl.n_tentative(:)';
n_c = fl.n_confirmed(:)';
n_d = fl.n_deleted(:)';
n_net = n_t + n_c - n_d;

out_dir = fullfile(this_dir, 'methods_v2_bayes', 'experiments', 'exp1_acoustic_ablation');
tex_dir = fullfile(this_dir, 'methods_v2_bayes', 'figures', 'ch3_exp1');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

fn = 'Microsoft YaHei';
fig = figure('Color','w','Position',[50 50 1040 640]);

ax1 = subplot(2,1,1);
hold on; grid on; box on;
plot(t, n_t, 'Color', [0.90 0.50 0.10], 'LineWidth', 1.6);
plot(t, n_c, 'Color', [0.15 0.40 0.80], 'LineWidth', 1.8);
plot(t, n_true, 'k-', 'LineWidth', 2.2);
ylabel('Count');
title(sprintf('Method B Field  seed %d: tentative / confirmed / GT', seed));
legend({'Tentative', 'Confirmed', 'GT'}, 'Location', 'northwest', 'Interpreter', 'none');
ylim([-0.3, max([n_t n_c n_true]) + 1.2]);
set(ax1, 'FontName', fn, 'FontSize', 11, 'XTickLabel', []);

ax2 = subplot(2,1,2);
hold on; grid on; box on;
plot(t, n_d, 'Color', [0.45 0.45 0.45], 'LineWidth', 1.8);
plot(t, n_net, 'Color', [0.55 0.15 0.55], 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Count');
title('Deleted (cumulative)  and  Tentative + Confirmed - Deleted');
legend({'Deleted', 'Tentative + Confirmed - Deleted'}, ...
    'Location', 'northwest', 'Interpreter', 'none');
set(ax2, 'FontName', fn, 'FontSize', 11);
linkaxes([ax1 ax2], 'x');
xlim([t(1), t(end)]);

exportgraphics(fig, fullfile(out_dir, 'fig_exp1_track_status.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(tex_dir, 'fig_exp1_track_status.png'), 'Resolution', 200);

fprintf('tentative  max=%d mean=%.2f\n', max(n_t), mean(n_t));
fprintf('confirmed  max=%d mean=%.2f\n', max(n_c), mean(n_c));
fprintf('deleted    end=%d\n', n_d(end));
fprintf('T+C-D      min=%d max=%d end=%d\n', min(n_net), max(n_net), n_net(end));
fprintf('saved fig_exp1_track_status.png\n');
