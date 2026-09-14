%% 同种子、无声学：只切换 pf_update_mode，验证 Radar-only 是否必然不同
clear; clc;

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

load(fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat'), 'allDetections', 'time');
data1 = load(fullfile(this_dir, 'Scenario', 'Step1_Data.mat'));
truth = data1.truth;

cfg0 = config_tracker();
cfg0.use_pf = true;
cfg0.save_mat = false;
cfg0.jpda_radar_only = true;
cfg0.birth_radar_only = true;
cfg0.birth_acoustic_eta = 0;

hs = [];
for ti = 1:numel(truth)
    zcol = truth(ti).pos3D(:,3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
ac_cfg = acoustic_config_from_scenario(data1.acousticPos, median(hs));
opts = struct('enable_acoustic', false, 'acoustic_mode', 'none', ...
    'verbose', false, 'method', 'softmax');

modes = {'legacy_likelihood', 'jpda_posterior'};
seeds = [1001, 42];
fprintf('seed | mode                 | OSPA   | FA%%   | card | peak | frag | medlife\n');
for si = 1:numel(seeds)
    for mi = 1:numel(modes)
        cfg = cfg0;
        cfg.pf_update_mode = modes{mi};
        rng(seeds(si), 'twister');
        out = run_tracker_once(allDetections, truth, time, cfg, ac_cfg, opts);
        m = eval_exp1_metrics(out.tracks, truth, time, out.cfg, out.frame_log, true);
        fprintf('%4d | %-20s | %6.2f | %5.1f | %4.2f | %4.0f | %4.0f | %6.1f\n', ...
            seeds(si), modes{mi}, m.ospa_mean, 100*m.false_confirm_ratio, ...
            m.card_err_mean, m.peak_confirmed, m.fragmentation_count, ...
            m.median_track_life);
    end
end
