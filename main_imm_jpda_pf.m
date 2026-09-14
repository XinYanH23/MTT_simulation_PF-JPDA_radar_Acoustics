%% =========================================================================
%  main_imm_jpda_pf.m
%  IMM-JPDA-PF 多目标融合跟踪主循环（声学概率场 + 粒子滤波）
%
%  与 main_imm_jpda_ekf.m 相比，本脚本仅修改：
%    1. cfg.use_pf = true（PF 替换 KF）
%    2. 每帧生成声学概率场 ac_pack_k
%    3. 历史数据记录（ac_pack_hist, assoc_hist, ess_log 用于可视化）
%    4. 调用 plot_results_pf 替代 plot_results
%
%  JPDA / IMM 模型概率 / 轨迹管理逻辑完全不变。
%% =========================================================================
clear; clc; close all;

%% ── 路径 ─────────────────────────────────────────────────────────────────
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

%% ── 参数 ──────────────────────────────────────────────────────────────────
cfg        = config_tracker();
cfg.use_pf = true;          % 开启 PF 模式（KF 所有模块自动回退，不影响 JPDA）

%% ── 数据加载 ──────────────────────────────────────────────────────────────
fprintf('===== 加载数据 =====\n');
data_dir_v3  = fullfile(this_dir, '..', 'v3');
data_dir_v2  = fullfile(this_dir, '..', 'v2.0');
scenario_dir = fullfile(this_dir, 'Scenario');

s2_cands = { fullfile(scenario_dir, 'Step2_HeteroDetections.mat'), ...
             fullfile(data_dir_v3,   'Step2_HeteroDetections_v3.mat'), ...
             fullfile(data_dir_v2,   'Step2_HeteroDetections.mat'),    ...
             'Step2_HeteroDetections.mat' };
s1_cands = { fullfile(scenario_dir, 'Step1_Data.mat'), ...
             fullfile(data_dir_v3,   'Step1_Data_v3.mat'), ...
             fullfile(data_dir_v2,   'Step1_Data.mat'),     ...
             'Step1_Data.mat' };

allDetections = [];  time = [];
for fi = 1:length(s2_cands)
    if exist(s2_cands{fi}, 'file')
        fprintf('  探测: %s\n', s2_cands{fi});
        load(s2_cands{fi}, 'allDetections', 'time');
        break
    end
end
data1 = [];
for fi = 1:length(s1_cands)
    if exist(s1_cands{fi}, 'file')
        fprintf('  场景: %s\n', s1_cands{fi});
        data1 = load(s1_cands{fi});
        break
    end
end
if isempty(allDetections) || isempty(data1)
    error('未找到数据文件，请先运行 run_scenario.m / step2_detection_simulation.m');
end

truth  = data1.truth;
N      = length(time);
cfg.dt = time(2) - time(1);
fprintf('  帧数: %d  目标: %d  dt: %.3f s\n', N, length(truth), cfg.dt);

%% ── 声学概率场配置（与场景 33 站对齐 + Step2 传播参数）────────────────
hs = [];
for ti = 1:length(truth)
    zcol = truth(ti).pos3D(:, 3);
    hs = [hs; zcol(~isnan(zcol))]; %#ok<AGROW>
end
if isempty(hs); src_h = 90; else; src_h = median(hs); end
if isfield(data1, 'acousticPos') && ~isempty(data1.acousticPos)
    ac_cfg = acoustic_config_from_scenario(data1.acousticPos, src_h);
    fprintf('  声学场: 场景 %d 站, src_h=%.1f m, A=%.0f, n=%.0f, nf=%.0f\n', ...
        size(data1.acousticPos,1), src_h, ac_cfg.A(1), ac_cfg.n(1), ac_cfg.noise_floor_db(1));
else
    ac_cfg = acoustic_config();
    fprintf('  声学场: 回退 acoustic_config() 默认 6 站\n');
end
prev_post_ac = [];

%% ── 初始化 ────────────────────────────────────────────────────────────────
tracks   = struct([]);
next_id  = 1;

ac_pack_hist = cell(N, 1);
assoc_hist   = cell(N, 1);

frame_log = struct('n_confirmed', zeros(1,N), ...
                   'n_tentative', zeros(1,N), ...
                   'n_deleted',   zeros(1,N));

snap_frames = round([N*0.25, N*0.5, N*0.75]);

fprintf('===== 开始 IMM-JPDA-PF 主循环 =====\n');

%% =========================================================================
%%  主循环
%% =========================================================================
for k = 1:N

    if mod(k, 50) == 0
        n_act = sum(~strcmp({tracks.status}, 'deleted'));
        fprintf('  帧 %4d/%d | 活跃: %d\n', k, N, n_act);
    end

    %% ── Step 0：声学概率场生成 ──────────────────────────────────────────
    ac_pack_k = [];
    src_pos_k = [];
    for tn = 1:length(truth)
        if k <= size(truth(tn).pos3D, 1) && ~any(isnan(truth(tn).pos3D(k,:)))
            src_pos_k = [src_pos_k; truth(tn).pos3D(k, 1:3)]; %#ok<AGROW>
        end
    end
    if ~isempty(src_pos_k)
        [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
            src_pos_k, ac_cfg, prev_post_ac, 'softmax');
        ac_pack_hist{k} = ac_pack_k;
    end

    %% ── Step 1：IMM 粒子混合 + 预测 + 融合 ─────────────────────────────
    for i = 1:length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end

        [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
                                                   tracks(i).mu, tracks(i).Pi);
        tracks(i).c_bar_ = c_bar_i;
        tracks(i).models = imm_pf_predict(tracks(i).models);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% ── Step 2&3：量测 + JPDA ──────────────────────────────────────────
    [Z_k, R_list_k, sensor_types_k] = detections_to_ZR(allDetections{k}, cfg);

    % 异构传感分离：JPDA 关联、PF 雷达似然与轨迹诞生仅使用雷达点迹；
    % 声学信息只经声学概率场进入 PF 权重（见 Step 4 的 ac_pack_k），避免重复计入。
    if (~isfield(cfg,'jpda_radar_only')) || cfg.jpda_radar_only
        if ~isempty(sensor_types_k)
            is_radar       = cellfun(@(s) strcmpi(s,'Radar'), sensor_types_k);
            Z_k            = Z_k(:, is_radar);
            R_list_k       = R_list_k(:, :, is_radar);
            sensor_types_k = sensor_types_k(is_radar);
        end
    end

    if isempty(tracks)
        active_mask = false(1, 0);
    else
        active_mask = ~strcmp({tracks.status}, 'deleted');
    end
    active_idx = find(active_mask);
    n_active   = length(active_idx);

    if n_active > 0 && size(Z_k, 2) > 0
        active_tracks = tracks(active_idx);
        assoc = jpda_run(active_tracks, Z_k, R_list_k, cfg);
    else
        assoc.valid_mat  = false(n_active, size(Z_k,2));
        assoc.innov_data = cell(n_active, size(Z_k,2));
        assoc.beta       = zeros(n_active, size(Z_k,2));
        assoc.beta0      = ones(n_active, 1);
        assoc.clusters   = {};
    end
    assoc_hist{k} = assoc;

    %% ── Step 4：PF 权重更新 ─────────────────────────────────────────────
    use_rp = isfield(cfg,'acoustic_assist_only') && cfg.acoustic_assist_only;
    for ti = 1:n_active
        i = active_idx(ti);
        if use_rp && ~isempty(sensor_types_k)
            is_radar = cellfun(@(s) strcmpi(s,'Radar'), sensor_types_k);
            if any(assoc.valid_mat(ti,:) & is_radar)
                assoc_ti = assoc;
                assoc_ti.valid_mat(ti, ~is_radar) = false;
                assoc_ti.beta(ti, ~is_radar)      = 0;
                assoc_ti.beta0(ti) = max(1-sum(assoc_ti.beta(ti,:)), 0);
            else
                assoc_ti = assoc;
            end
        else
            assoc_ti = assoc;
        end
        tracks(i) = track_imm_pda_update(tracks(i), Z_k, R_list_k, ...
                                          assoc_ti, ti, cfg, ac_pack_k);

        if ismember(k, snap_frames)
            si = find(snap_frames == k, 1);
            [~, mj] = max(tracks(i).mu);
            tracks(i).particle_snap{si} = tracks(i).models(mj).particles(1:2, :);
        end
    end

    %% ── Step 5：轨迹管理（诞生过滤 + 并发冗余去重）──────────────────────
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, ...
                                      k, next_id, cfg, sensor_types_k, ac_pack_k);
    tracks = track_merge(tracks, cfg);

    %% ── 帧统计 ────────────────────────────────────────────────────────
    if ~isempty(tracks)
        st = {tracks.status};
        frame_log.n_confirmed(k) = sum(strcmp(st,'confirmed'));
        frame_log.n_tentative(k) = sum(strcmp(st,'tentative'));
        frame_log.n_deleted(k)   = sum(strcmp(st,'deleted'));
    end
    if ~isempty(tracks) && isfield(tracks,'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end
end

fprintf('===== 主循环完成 =====\n');

%% ── 性能评估 ──────────────────────────────────────────────────────────────
metrics = eval_metrics(tracks, truth, time, cfg);

%% ── 保存 ──────────────────────────────────────────────────────────────────
if cfg.save_mat
    out_path = fullfile(this_dir, 'Step3_IMM_JPDA_PF_Result.mat');
    save(out_path, 'tracks','metrics','frame_log','truth','time','cfg', ...
         'ac_pack_hist','assoc_hist', '-v7.3');
    fprintf('结果已保存: %s\n', out_path);
end

%% ── 可视化（8 幅图） ──────────────────────────────────────────────────────
fig_pf = plot_results_pf(tracks, truth, time, metrics, cfg, ac_pack_hist, assoc_hist);
fig_path = fullfile(this_dir, 'report', 'IMM_JPDA_PF_Results.png');
if ~exist(fullfile(this_dir, 'report'), 'dir')
    mkdir(fullfile(this_dir, 'report'));
end
if ~isempty(fig_pf) && ishandle(fig_pf)
    saveas(fig_pf, fig_path);
    fprintf('可视化已保存: %s\n', fig_path);
end
fprintf('===== 完成 =====\n');
