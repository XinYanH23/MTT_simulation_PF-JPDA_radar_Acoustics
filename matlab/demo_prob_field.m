function demo_prob_field()
%DEMO_PROB_FIELD 声学概率场模块端到端演示（无需音频）。
%
%   流程：
%     1. 单节点探测能力随距离曲线 + 地面反射敏感性
%     2. 仿真一条声源轨迹的各节点 SPL
%     3. 逐帧构造概率场（三种映射）+ 递归贝叶斯滤波
%     4. 统计 MAP/均值/瞬时 三种点估计误差并作图
%
%   直接运行：>> demo_prob_field

cfg = acoustic_config();
fprintf('=== 声学概率场模块演示 ===\n');

% ---------------------------------------------------------------------
% 1. 探测能力随距离 + 地面反射敏感性
% ---------------------------------------------------------------------
figure('Name', '探测能力随距离', 'Color', 'w');
for m = 1:size(cfg.mic_coords,1)
    cap = detection_capability(m, cfg);
    subplot(1,2,1); plot(cap.r_m, cap.snr_db, 'LineWidth', 1.4); hold on;
    subplot(1,2,2); plot(cap.r_m, cap.p_detect, 'LineWidth', 1.4); hold on;
    fprintf('mic%d 探测半径 r(P_D=0.5) = %.1f m\n', m, cap.r50_m);
end
subplot(1,2,1); grid on; xlabel('距离 r (m)'); ylabel('SNR (dB)'); title('SNR-距离'); legend(cfg.mic_ids);
subplot(1,2,2); grid on; xlabel('距离 r (m)'); ylabel('P_D'); title('探测概率-距离'); legend(cfg.mic_ids);

% 地面反射敏感性（mic1）
figure('Name', '地面反射敏感性', 'Color', 'w');
Rs = [0.0, 0.3, 0.6, 0.9];
for R = Rs
    c2 = cfg; c2.R_ground = R;
    cap = detection_capability(1, c2);
    plot(cap.r_m, cap.Lp_db, 'LineWidth', 1.4); hold on;
end
grid on; xlabel('距离 r (m)'); ylabel('Lp (dB)');
title('mic1 不同地面反射系数 R 下的 SPL-距离'); legend(arrayfun(@(r) sprintf('R=%.1f',r), Rs, 'uni', 0));

% ---------------------------------------------------------------------
% 2. 仿真声源轨迹
% ---------------------------------------------------------------------
rng(42);
T = 40;
t = linspace(0, 1, T);
% 在麦克风包络内画一条弧线轨迹
cx = mean(cfg.mic_coords(:,1)); cz = mean(cfg.mic_coords(:,3));
traj = [cx + 6000*cos(2*pi*t') - 3000, ...
        repmat(cfg.src_height_mm, T, 1), ...
        cz + 6000*sin(2*pi*t')];

% ---------------------------------------------------------------------
% 3. 逐帧概率场 + 递归贝叶斯
% ---------------------------------------------------------------------
err_map = nan(T,1); err_mean = nan(T,1); err_inst = nan(T,1);
prev_post = [];
for k = 1:T
    measured_lp = simulate_node_spl(traj(k,:), cfg, cfg.sigma0_db);
    pack = build_acoustic_prob_field(measured_lp, cfg, 'softmax', prev_post);
    prev_post = pack.post;

    gt = [traj(k,1), traj(k,3)];
    err_map(k)  = norm(pack.est_map  - gt) / 1000;   % m
    err_mean(k) = norm(pack.est_mean - gt) / 1000;
    % 瞬时（无时间约束）：直接对似然场取均值
    [im, ~] = field_point_estimate(pack.L_softmax, pack.XX, pack.ZZ);
    err_inst(k) = norm(im - gt) / 1000;
end

fprintf('\n定位误差 (m)：瞬时均值=%.2f  贝叶斯MAP=%.2f  贝叶斯均值=%.2f\n', ...
    nanmean(err_inst), nanmean(err_map), nanmean(err_mean));

% ---------------------------------------------------------------------
% 4. 可视化最后一帧概率场 + 误差曲线
% ---------------------------------------------------------------------
figure('Name', '概率场（末帧）', 'Color', 'w');
fields = {pack.L_softmax, pack.L_intensity, pack.L_nms_gmm, pack.post};
names  = {'方案一 Softmax', '方案二 能量域', '方案三 NMS-GMM', '递归贝叶斯后验'};
for i = 1:4
    subplot(2,2,i);
    imagesc(pack.zc/1000, pack.xc/1000, fields{i}); axis xy; hold on;
    plot(traj(end,3)/1000, traj(end,1)/1000, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(cfg.mic_coords(:,3)/1000, cfg.mic_coords(:,1)/1000, 'w^', 'MarkerFaceColor','k');
    xlabel('Z (m)'); ylabel('X (m)'); title(names{i}); colorbar;
end

figure('Name', '逐帧定位误差', 'Color', 'w');
plot(err_inst, 'r--', 'LineWidth', 1.2); hold on;
plot(err_map, 'm-', 'LineWidth', 1.2);
plot(err_mean, 'g-', 'LineWidth', 1.6);
grid on; xlabel('帧'); ylabel('误差 (m)');
legend('瞬时均值', '贝叶斯 MAP', '贝叶斯均值'); title('三种估计策略逐帧误差');

fprintf('演示完成。\n');
end
