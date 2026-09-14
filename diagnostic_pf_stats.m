%% diagnostic_pf_stats.m
%  PF 粒子健康度诊断脚本
%  不修改任何现有代码；直接调用已有函数复现主循环前50帧，采集粒子统计量。
%
%  输出：
%    1. 每帧 std(px), std(py)（前50帧，所有活跃轨迹取均值）
%    2. ESS 历史：min / mean / max
%    3. 重采样次数 resample_count
%    4. 预测步 diag(Q) 两个模型
%    5. roughening/jitter 是否存在（静态结论）
%    6. 唯一粒子比例 unique_particles/N（前50帧）
%    7. 折线图：粒子位置标准差随帧变化

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

%% ── 参数 ──────────────────────────────────────────────────────────────────
cfg        = config_tracker();
cfg.use_pf = true;
ac_cfg     = acoustic_config();
prev_post_ac = [];

%% ── 数据加载 ──────────────────────────────────────────────────────────────
scenario_dir = fullfile(this_dir, 'Scenario');
s2 = fullfile(scenario_dir, 'Step2_HeteroDetections.mat');
s1 = fullfile(scenario_dir, 'Step1_Data.mat');
if ~exist(s2,'file') || ~exist(s1,'file')
    error('数据文件未找到，请先运行 run_scenario.m / step2_detection_simulation.m');
end
load(s2, 'allDetections', 'time');
data1 = load(s1);
truth = data1.truth;
N_total = length(time);
cfg.dt  = time(2) - time(1);

N_DIAG  = min(50, N_total);   % 只跑前50帧

fprintf('===== PF 诊断脚本（前 %d 帧） =====\n', N_DIAG);
fprintf('  dt = %.4f s,  pf_N = %d,  resample_thresh = %.2f\n', ...
    cfg.dt, cfg.pf_N, cfg.pf_resample_thresh);

%% ── 第4项：Q 对角元素（静态计算）─────────────────────────────────────────
dt = cfg.dt;
for jm = 1:cfg.n_models
    qv = [cfg.cv_q, cfg.ca_q]; q = qv(jm);
    Q_diag = [q*dt^3/3, q*dt^3/3, q*dt, q*dt];
    fprintf('\n[Q model%d  q=%.1f  dt=%.4f]\n', jm, q, dt);
    fprintf('  diag(Q) = [%.4e, %.4e, %.4e, %.4e]\n', Q_diag);
    fprintf('  sqrt(diag(Q)) = [%.4f m, %.4f m, %.4f m/s, %.4f m/s]\n', sqrt(Q_diag));
end

%% ── 第5项：roughening/jitter 检查（静态）─────────────────────────────────
fprintf('\n[Roughening/Jitter]\n');
fprintf('  pf_systematic_resample.m 中仅做 idx 重索引，无噪声注入。\n');
fprintf('  工程中 *无* 显式 roughening。\n');
fprintf('  隐式扩散：imm_pf_predict.m 中 L*randn(4,N) 在重采样后下一帧仍注入过程噪声。\n');

%% ── 初始化主循环状态 ──────────────────────────────────────────────────────
tracks    = struct([]);
next_id   = 1;

% 诊断记录（每帧）
diag_px_std      = nan(N_DIAG, 1);
diag_py_std      = nan(N_DIAG, 1);
diag_ess         = nan(N_DIAG, 1);
diag_unique_ratio= nan(N_DIAG, 1);
resample_count   = 0;

%% ── 包装 pf_systematic_resample 以计数（闭包变量）──────────────────────
% 使用全局变量计数（不修改任何原有函数）
global DIAG_RESAMPLE_COUNT
DIAG_RESAMPLE_COUNT = 0;

%% ── 主循环（前 N_DIAG 帧）────────────────────────────────────────────────
snap_frames = round([N_total*0.25, N_total*0.5, N_total*0.75]);

for k = 1:N_DIAG

    %% Step 0：声学概率场
    ac_pack_k = [];
    src_pos_k = [];
    for tn = 1:length(truth)
        if k <= size(truth(tn).pos3D, 1) && ~any(isnan(truth(tn).pos3D(k,:)))
            src_pos_k = [src_pos_k; truth(tn).pos3D(k, 1:3)]; %#ok
        end
    end
    if ~isempty(src_pos_k)
        try
            [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
                src_pos_k, ac_cfg, prev_post_ac, 'softmax');
        catch ME
            warning('acoustic_field_to_tracker 失败帧%d: %s', k, ME.message);
        end
    end

    %% Step 1：IMM 预测
    for i = 1:length(tracks)
        if strcmp(tracks(i).status, 'deleted'), continue; end
        [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
                                                   tracks(i).mu, tracks(i).Pi);
        tracks(i).c_bar_ = c_bar_i;
        tracks(i).models = imm_pf_predict(tracks(i).models);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% Step 2&3：量测 + JPDA
    [Z_k, R_list_k, ~] = detections_to_ZR(allDetections{k}, cfg);

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

    %% Step 4：PF 权重更新
    for ti = 1:n_active
        i = active_idx(ti);
        tracks(i) = track_imm_pda_update(tracks(i), Z_k, R_list_k, ...
                                          assoc, ti, cfg, ac_pack_k);
        if ismember(k, snap_frames)
            si = find(snap_frames == k, 1);
            [~, mj] = max(tracks(i).mu);
            tracks(i).particle_snap{si} = tracks(i).models(mj).particles(1:2,:);
        end
    end

    %% Step 5：轨迹管理
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, ...
                                      k, next_id, cfg);

    if ~isempty(tracks) && isfield(tracks, 'c_bar_')
        tracks = rmfield(tracks, 'c_bar_');
    end

    %% ── 诊断统计采集 ──────────────────────────────────────────────────────
    if ~isempty(tracks)
        all_px = [];  all_py = [];
        ess_vals = [];  n_unique_total = 0;  n_total = 0;

        active_now = find(~strcmp({tracks.status}, 'deleted'));
        for ii = active_now
            for jm = 1:cfg.n_models
                pts = tracks(ii).models(jm).particles;  % 4×N
                w   = tracks(ii).models(jm).weights;    % 1×N

                % 粒子位置标准差（粒子云离散度，非加权，反映多样性）
                all_px = [all_px, pts(1,:)]; %#ok
                all_py = [all_py, pts(2,:)]; %#ok

                % ESS
                ess_now = 1 / sum(w.^2);
                ess_vals = [ess_vals, ess_now]; %#ok

                % 唯一粒子比例（取前两维 px,py 判重）
                pts_int = round(pts(1:2,:) * 10);  % 精度 0.1m
                n_unique = size(unique(pts_int', 'rows'), 1);
                n_unique_total = n_unique_total + n_unique;
                n_total        = n_total + cfg.pf_N;

                % 重采样计数（检测是否本帧发生过重采样：ESS < thresh 则触发）
                if ess_now < cfg.pf_resample_thresh * cfg.pf_N
                    resample_count = resample_count + 1;
                end
            end
        end

        if ~isempty(all_px)
            diag_px_std(k) = std(all_px);
            diag_py_std(k) = std(all_py);
        end
        if ~isempty(ess_vals)
            diag_ess(k) = mean(ess_vals);
        end
        if n_total > 0
            diag_unique_ratio(k) = n_unique_total / n_total;
        end
    end

end  % 主循环

%% ── 输出第1项：逐帧 std(px), std(py) ────────────────────────────────────
fprintf('\n========== 第1项：逐帧粒子位置标准差（前%d帧）==========\n', N_DIAG);
fprintf('%-6s  %-12s  %-12s\n', '帧', 'std(px)/m', 'std(py)/m');
fprintf('%s\n', repmat('-',1,35));
for k = 1:N_DIAG
    if ~isnan(diag_px_std(k))
        fprintf('  %3d    %10.4f    %10.4f\n', k, diag_px_std(k), diag_py_std(k));
    end
end

%% ── 输出第2项：ESS 统计 ──────────────────────────────────────────────────
valid_ess = diag_ess(~isnan(diag_ess));
fprintf('\n========== 第2项：ESS 统计 ==========\n');
if ~isempty(valid_ess)
    fprintf('  ESS min  = %.2f  (N=%d, 占比 %.1f%%)\n', min(valid_ess),  cfg.pf_N, min(valid_ess)/cfg.pf_N*100);
    fprintf('  ESS mean = %.2f  (N=%d, 占比 %.1f%%)\n', mean(valid_ess), cfg.pf_N, mean(valid_ess)/cfg.pf_N*100);
    fprintf('  ESS max  = %.2f  (N=%d, 占比 %.1f%%)\n', max(valid_ess),  cfg.pf_N, max(valid_ess)/cfg.pf_N*100);
    thresh_frames = sum(valid_ess < cfg.pf_resample_thresh * cfg.pf_N);
    fprintf('  ESS < 阈值(%.0f)的帧数: %d / %d (%.1f%%)\n', ...
        cfg.pf_resample_thresh*cfg.pf_N, thresh_frames, length(valid_ess), ...
        thresh_frames/length(valid_ess)*100);
else
    fprintf('  无有效 ESS 数据（前%d帧内无活跃轨迹）\n', N_DIAG);
end

%% ── 输出第3项：重采样次数 ────────────────────────────────────────────────
fprintf('\n========== 第3项：重采样次数 ==========\n');
fprintf('  resample_count（ESS<阈值触发次数）= %d\n', resample_count);
fprintf('  [注] 此计数为 ESS检测次数，与 pf_systematic_resample 实际调用次数等价\n');
fprintf('       (pf_jpda_update.m:120 if ess < thresh → pf_systematic_resample)\n');

%% ── 输出第6项：唯一粒子比例 ──────────────────────────────────────────────
fprintf('\n========== 第6项：唯一粒子比例（前%d帧）==========\n', N_DIAG);
fprintf('%-6s  %-20s\n', '帧', 'unique_particles/N');
fprintf('%s\n', repmat('-',1,30));
for k = 1:N_DIAG
    if ~isnan(diag_unique_ratio(k))
        fprintf('  %3d    %6.3f\n', k, diag_unique_ratio(k));
    end
end
valid_ur = diag_unique_ratio(~isnan(diag_unique_ratio));
if ~isempty(valid_ur)
    fprintf('\n  唯一粒子比例：min=%.3f  mean=%.3f  max=%.3f\n', ...
        min(valid_ur), mean(valid_ur), max(valid_ur));
end

%% ── 第7项：绘图 ──────────────────────────────────────────────────────────
frames_with_data = find(~isnan(diag_px_std));

fig = figure('Name', 'PF 粒子诊断', 'Position', [100 100 1200 500], 'Color', 'w');

%% 子图1：粒子位置标准差
ax1 = subplot(1, 3, 1); hold on; grid on; box on;
title('粒子位置标准差', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('帧'); ylabel('std (m)');
if ~isempty(frames_with_data)
    plot(frames_with_data, diag_px_std(frames_with_data), 'b-o', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'std(px)');
    plot(frames_with_data, diag_py_std(frames_with_data), 'r-s', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'std(py)');
    legend('Location', 'best', 'FontSize', 9);
    ylim_max = max([diag_px_std(frames_with_data); diag_py_std(frames_with_data)]) * 1.15;
    if isfinite(ylim_max) && ylim_max > 0
        ylim([0, ylim_max]);
    end
else
    text(0.5, 0.5, '无数据', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlim([1, N_DIAG]);

%% 子图2：ESS 变化曲线
ax2 = subplot(1, 3, 2); hold on; grid on; box on;
title('ESS 随帧变化', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('帧'); ylabel('ESS');
frames_ess = find(~isnan(diag_ess));
if ~isempty(frames_ess)
    plot(frames_ess, diag_ess(frames_ess), 'm-o', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'ESS(mean)');
    yline(cfg.pf_resample_thresh * cfg.pf_N, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('阈值 %.0f', cfg.pf_resample_thresh*cfg.pf_N));
    yline(cfg.pf_N, 'k:', 'LineWidth', 1, 'DisplayName', sprintf('N=%d', cfg.pf_N));
    legend('Location', 'best', 'FontSize', 9);
    ylim([0, cfg.pf_N * 1.1]);
else
    text(0.5, 0.5, '无数据', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlim([1, N_DIAG]);

%% 子图3：唯一粒子比例
ax3 = subplot(1, 3, 3); hold on; grid on; box on;
title('唯一粒子比例  unique/N', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('帧'); ylabel('unique\_particles / N');
frames_ur = find(~isnan(diag_unique_ratio));
if ~isempty(frames_ur)
    plot(frames_ur, diag_unique_ratio(frames_ur), 'g-^', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'unique ratio');
    yline(1.0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    legend('Location', 'best', 'FontSize', 9);
    ylim([0, 1.05]);
else
    text(0.5, 0.5, '无数据', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlim([1, N_DIAG]);

sgtitle(sprintf('PF 粒子诊断  |  N=%d  |  前%d帧', cfg.pf_N, N_DIAG), ...
    'FontSize', 12, 'FontWeight', 'bold');
drawnow;

% 保存图像
fig_path = fullfile(this_dir, 'report', 'pf_diagnostic.png');
if ~exist(fullfile(this_dir, 'report'), 'dir'); mkdir(fullfile(this_dir, 'report')); end
saveas(fig, fig_path);
fprintf('\n诊断图已保存: %s\n', fig_path);
fprintf('===== 诊断完成 =====\n');
