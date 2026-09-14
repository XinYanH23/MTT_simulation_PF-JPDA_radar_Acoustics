%% EXPORT_EXP3_FOR_TEX  从蒙特卡洛 exp3_results.mat 重绘
clear; clc; close all;
exp_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(exp_dir));
addpath(fullfile(root, 'Evaluation'));
out_dir = exp_dir;
tex_dir = fullfile(root, 'report', 'figures', 'ch3_exp3');
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end

S = load(fullfile(exp_dir, 'exp3_results.mat'));
assert(isfield(S.metrics_all, 'mu0_RadarOnly') && isfield(S.metrics_all.mu0_RadarOnly, 'agg'), ...
    '请先运行新版 run_exp3_clutter_sweep.m（蒙特卡洛）');
fprintf('Exp3 export: 请直接运行 run_exp3_clutter_sweep.m，其已写出图与 tab_exp3_metrics.tex\n');
fprintf('mat: %s  n_mc=%d\n', fullfile(exp_dir, 'exp3_results.mat'), S.n_mc);
