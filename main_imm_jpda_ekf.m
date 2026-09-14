%% =========================================================================
%  main_imm_jpda_ekf.m
%  IMM-JPDA-EKF 多目标融合跟踪基线系统 — 主循环
%
%  架构（Bar-Shalom 2001 + Blackman & Popoli 1999 经典流程）：
%
%    Step1_Data.mat / Step2_HeteroDetections.mat
%          ↓
%    [每帧循环]
%      1. IMM 预测（混合 + 各模型 EKF 预测 + 融合）
%      2. 量测门控（马氏距离验证门）
%      3. JPDA 关联（联合事件枚举 + β 边缘化）
%      4. PDA 更新（完整三项式协方差，含 innovation spread）
%      5. 轨迹确认（M/N 逻辑）
%      6. 轨迹诞生（未关联探测 → 新 Tentative 轨迹）
%      7. 轨迹删除（连续漏检超限）
%      8. 轨迹输出
%          ↓
%    性能评估 + 可视化
%
%  使用方法：
%    1. 将 MATLAB 当前路径设置为 Tracker/ 目录（或添加到 path）
%    2. 确保 Step2_HeteroDetections.mat 和 Step1_Data.mat 可访问
%    3. 直接运行本脚本
% =========================================================================
clear; clc; close all;

%% ── 路径设置 ────────────────────────────────────────────────────────────
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

%% ── 参数配置 ────────────────────────────────────────────────────────────
cfg = config_tracker();

% 声学概率场配置（KF/PF 共用同一场输入；use_pf 只切换滤波器）
ac_cfg = acoustic_config();
prev_post_ac = [];

%% ── 数据加载 ────────────────────────────────────────────────────────────
fprintf('===== 加载数据 =====\n');

data_dir_v3  = fullfile(this_dir, '..', 'v3');
data_dir_v2  = fullfile(this_dir, '..', 'v2.0');
scenario_dir = fullfile(this_dir, 'Scenario');

% 优先 Scenario/，其次 v3，再回退 v2.0
s2_candidates = { fullfile(scenario_dir, 'Step2_HeteroDetections.mat'), ...
                  fullfile(data_dir_v3, 'Step2_HeteroDetections_v3.mat'), ...
                  fullfile(data_dir_v2, 'Step2_HeteroDetections.mat'),    ...
                  'Step2_HeteroDetections.mat' };
s1_candidates = { fullfile(scenario_dir, 'Step1_Data.mat'), ...
                  fullfile(data_dir_v3, 'Step1_Data_v3.mat'), ...
                  fullfile(data_dir_v2, 'Step1_Data.mat'),     ...
                  'Step1_Data.mat' };

allDetections = [];  time = [];
for fi = 1 : length(s2_candidates)
    if exist(s2_candidates{fi}, 'file')
        fprintf('  探测数据: %s\n', s2_candidates{fi});
        load(s2_candidates{fi}, 'allDetections', 'time');
        break
    end
end

data1 = [];
for fi = 1 : length(s1_candidates)
    if exist(s1_candidates{fi}, 'file')
        fprintf('  场景数据: %s\n', s1_candidates{fi});
        data1 = load(s1_candidates{fi});
        break
    end
end

if isempty(allDetections) || isempty(data1)
    error('未找到所需 .mat 数据文件，请先运行 Step1/Step2 生成数据。');
end

truth    = data1.truth;
N        = length(time);
cfg.dt   = time(2) - time(1);

fprintf('  帧数: %d,  真值目标数: %d,  采样间隔: %.3f s\n', N, length(truth), cfg.dt);

%% ── 初始化 ──────────────────────────────────────────────────────────────
tracks  = struct([]);   % 0×0 结构体数组（避免 [] 导致后续字段不一致）
next_id = 1;            % 全局轨迹 ID 计数器

% 预分配性能记录
frame_log = struct('n_confirmed', zeros(1,N), ...
                   'n_tentative', zeros(1,N), ...
                   'n_deleted',   zeros(1,N));

fprintf('===== 开始主跟踪循环 =====\n');

%% =========================================================================
%% 主循环
%% =========================================================================
for k = 1 : N

    if mod(k, 50) == 0
        n_act = sum(~strcmp({tracks.status}, 'deleted'));
        fprintf('  帧 %4d/%d | 活跃轨迹: %d\n', k, N, n_act);
    end

    %% ── 步骤 0b：声学概率场生成（每帧；KF/PF 共用） ─────────────────────
    ac_pack_k = [];
    if ~isempty(ac_cfg)
        % 提取当前帧真值位置（供声学仿真使用）
        src_pos_k = [];
        for tn = 1:length(truth)
            if k <= size(truth(tn).pos3D, 1) && ~any(isnan(truth(tn).pos3D(k,:)))
                src_pos_k = [src_pos_k; truth(tn).pos3D(k, 1:3)]; %#ok<AGROW>
            end
        end
        if ~isempty(src_pos_k)
            [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
                src_pos_k, ac_cfg, prev_post_ac, 'softmax');
        end
    end

    %% ── 步骤 1：IMM 预测（Mixing + Predict + Fuse） ─────────────────────
    for i = 1 : length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end

        if isfield(cfg, 'use_pf') && cfg.use_pf
            % 1a. PF 粒子混合
            [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
                                                       tracks(i).mu, ...
                                                       tracks(i).Pi);
            tracks(i).c_bar_ = c_bar_i;
            % 1b. 粒子预测
            tracks(i).models = imm_pf_predict(tracks(i).models);
        else
            % 1a. KF IMM 混合（原有，不变）
            [tracks(i).models, c_bar_i] = imm_mix(tracks(i).models, ...
                                                   tracks(i).mu, ...
                                                   tracks(i).Pi);
            tracks(i).c_bar_ = c_bar_i;
            % 1b. KF 预测（原有，不变）
            tracks(i).models = imm_predict(tracks(i).models);
        end

        % 1c. IMM 融合预测状态（供门控；KF/PF 均调用相同 imm_fuse）
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% ── 步骤 2 & 3：量测获取 + JPDA 关联 ───────────────────────────────
    [Z_k, R_list_k, sensor_types_k] = detections_to_ZR(allDetections{k}, cfg);

    % 提取活跃轨迹（tracks=[] 时不能对 status 做点索引）
    if isempty(tracks)
        active_mask = false(1, 0);
    else
        active_mask = ~strcmp({tracks.status}, 'deleted');
    end
    active_idx  = find(active_mask);
    n_active    = length(active_idx);

    if n_active > 0 && size(Z_k, 2) > 0
        active_tracks = tracks(active_idx);
        assoc = jpda_run(active_tracks, Z_k, R_list_k, cfg);
    else
        % 无活跃轨迹或无量测：构造空 assoc
        assoc.valid_mat  = false(n_active, size(Z_k,2));
        assoc.innov_data = cell(n_active, size(Z_k,2));
        assoc.beta       = zeros(n_active, size(Z_k,2));
        assoc.beta0      = ones(n_active, 1);
        assoc.clusters   = {};
    end

    %% ── 步骤 4：PDA 更新（每条活跃轨迹）────────────────────────────────
    % 雷达优先融合：若 cfg.acoustic_assist_only=true，对每条轨迹独立过滤
    % 门内量测集合，保留雷达探测；仅在无雷达探测时回退使用声学探测。
    % 此过滤仅作用于 PDA 更新，不修改 assoc（轨迹管理仍使用原始 assoc）。
    use_radar_priority = isfield(cfg, 'acoustic_assist_only') && cfg.acoustic_assist_only;

    for ti = 1 : n_active
        i = active_idx(ti);

        if use_radar_priority && ~isempty(sensor_types_k)
            % 找出本轨迹门内的雷达探测列索引
            is_radar = cellfun(@(s) strcmpi(s, 'Radar'), sensor_types_k);
            radar_in_gate = assoc.valid_mat(ti, :) & is_radar;

            if any(radar_in_gate)
                % 有雷达探测：屏蔽声学探测，重新归一化 beta
                assoc_ti          = assoc;
                acoustic_mask     = ~is_radar;
                assoc_ti.valid_mat(ti, acoustic_mask) = false;
                assoc_ti.beta(ti,  acoustic_mask)     = 0;
                % beta 之和已减少，相应提高 beta0 使概率守恒
                sum_beta_ti = sum(assoc_ti.beta(ti, :));
                assoc_ti.beta0(ti) = max(1 - sum_beta_ti, 0);
            else
                % 无雷达探测：声学回退，不改变 assoc
                assoc_ti = assoc;
            end
        else
            assoc_ti = assoc;
        end

        tracks(i) = track_imm_pda_update(tracks(i), Z_k, R_list_k, ...
                                          assoc_ti, ti, cfg, ac_pack_k);
    end

    %% ── 步骤 5~7：轨迹管理（确认 / 诞生 / 删除）────────────────────────
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, ...
                                      k, next_id, cfg);

    %% ── 步骤 8：帧统计记录 ──────────────────────────────────────────────
    if ~isempty(tracks)
        statuses = {tracks.status};
        frame_log.n_confirmed(k) = sum(strcmp(statuses, 'confirmed'));
        frame_log.n_tentative(k) = sum(strcmp(statuses, 'tentative'));
        frame_log.n_deleted(k)   = sum(strcmp(statuses, 'deleted'));
    end

    % 清理临时字段（须整数组 rmfield，避免结构体字段不一致）
    if ~isempty(tracks) && isfield(tracks, 'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end
end   % 主循环结束

fprintf('===== 跟踪循环完成 =====\n');

%% ── 性能评估 ────────────────────────────────────────────────────────────
metrics = eval_metrics(tracks, truth, time, cfg);

%% ── 结果保存 ────────────────────────────────────────────────────────────
if cfg.save_mat
    out_path = fullfile(this_dir, 'Step3_IMM_JPDA_EKF_Result.mat');
    save(out_path, 'tracks', 'metrics', 'frame_log', 'truth', 'time', 'cfg', '-v7.3');
    fprintf('结果已保存: %s\n', out_path);
end

%% ── 可视化 ──────────────────────────────────────────────────────────────
if isfield(cfg, 'use_pf') && cfg.use_pf
    plot_results_pf(tracks, truth, time, metrics, cfg);
else
    plot_results(tracks, truth, time, metrics, cfg);
end
if exist(fullfile(this_dir, 'Scenario', 'plot_full_pipeline.m'), 'file')
    try
        plot_full_pipeline();
    catch ME
        warning('plot_full_pipeline:Failed', '%s', ME.message);
    end
end

fprintf('===== 完成 =====\n');
