function demo_mic_placement()
%DEMO_MIC_PLACEMENT 贪婪麦克风布设算法端到端演示。
%
%   演示：在 20 个待选点中，用贪婪算法选出满足全域平均 P_D_joint >= 0.80 的
%   最少节点布局，并可视化每步边际增益与最终覆盖图。

cfg = acoustic_config();
rng(7);

% ------------------------------------------------------------------
% 1. 生成待选点（示例：在操作区域内均匀散布 20 个候选地面节点）
% ------------------------------------------------------------------
% 覆盖评估区域（与三个真实麦克风包络对齐）
pts = cfg.mic_coords;
xmin = min(pts(:,1)) - 5000; xmax = max(pts(:,1)) + 5000;
zmin = min(pts(:,3)) - 5000; zmax = max(pts(:,3)) + 5000;
cfg.placement_area_mm = [xmin, xmax, zmin, zmax];

K = 20;
cx = xmin + rand(K,1)*(xmax-xmin);
cz = zmin + rand(K,1)*(zmax-zmin);
cy = zeros(K,1);                    % 地面高度 Y=0
candidates = [cx, cy, cz];

fprintf('待选点数 K=%d，覆盖区域 [%.0f~%.0f x %.0f~%.0f] mm\n', ...
    K, xmin, xmax, zmin, zmax);

% ------------------------------------------------------------------
% 2. 运行贪婪布设算法
% ------------------------------------------------------------------
result = greedy_mic_placement(candidates, cfg);

fprintf('\n贪婪布设结果：\n');
fprintf('  选取节点数  : %d\n', result.n_nodes);
fprintf('  达到覆盖目标: %d (目标 mean_pd >= %.2f)\n', result.satisfied, cfg.placement_min_coverage);
fprintf('  最终 mean_pd : %.4f\n', result.mean_pd_hist(end));
fprintf('  最终 coverage_frac (P_D>=0.5): %.4f\n', result.coverage_hist(end));
fprintf('  选取顺序（候选点行号）: ');
fprintf('%d ', result.sel_idx); fprintf('\n');
fprintf('  边际增益：'); fprintf('%.4f ', result.marginal_gain); fprintf('\n');

% ------------------------------------------------------------------
% 3. 对比：空集 / 逐步 / 全部 20 个节点的覆盖
% ------------------------------------------------------------------
[pd_all, xc, zc] = compute_coverage_map(candidates, ...
    repmat(cfg.placement_A_default, K,1), ...
    repmat(cfg.placement_n_default, K,1), ...
    repmat(cfg.placement_nf_default, K,1), cfg);
fprintf('  使用全部 %d 个候选节点的 mean_pd = %.4f\n', K, mean(pd_all(:)));

% ------------------------------------------------------------------
% 4. 可视化
% ------------------------------------------------------------------
figure('Name','贪婪布设：逐步覆盖率', 'Color','w');
steps = 0:result.n_nodes;
yyaxis left;
plot(steps, result.mean_pd_hist, 'bo-', 'LineWidth', 1.6);
yline(cfg.placement_min_coverage, 'b--', sprintf('目标 %.2f', cfg.placement_min_coverage));
ylabel('全域平均 P_D\_joint');
yyaxis right;
bar(1:result.n_nodes, result.marginal_gain, 0.4, 'FaceColor', [0.8 0.4 0.2], 'FaceAlpha', 0.6);
ylabel('边际增益');
xlabel('已选节点数'); grid on;
title('贪婪布设：逐步全域平均 P_D 与边际增益');

figure('Name','最终覆盖图', 'Color','w');
imagesc(result.zc/1000, result.xc/1000, result.pd_map_final);
axis xy; hold on; colorbar; clim([0 1]);
xlabel('Z (m)'); ylabel('X (m)'); title('最终 P\_D\_joint 分布');
% 画出所有候选（灰圆）与已选（红星）
plot(candidates(:,3)/1000, candidates(:,1)/1000, 'o', ...
    'MarkerSize',6, 'Color',[0.6 0.6 0.6], 'MarkerFaceColor',[0.8 0.8 0.8]);
sel = result.sel_coords;
for i = 1:result.n_nodes
    plot(sel(i,3)/1000, sel(i,1)/1000, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    text(sel(i,3)/1000+0.2, sel(i,1)/1000, sprintf('#%d',i), 'Color','r', 'FontSize',9);
end
% 画出真实麦克风位置（白三角）
plot(cfg.mic_coords(:,3)/1000, cfg.mic_coords(:,1)/1000, 'w^', ...
    'MarkerFaceColor','k', 'MarkerSize',10);
legend('候选点','已选节点','真实麦克风', 'Location','best');

% ------------------------------------------------------------------
% 5. 与真实三麦克风的覆盖对比
% ------------------------------------------------------------------
[pd_real, ~, ~] = compute_coverage_map(cfg.mic_coords, cfg.A, cfg.n, cfg.noise_floor_db, cfg);
fprintf('\n对比（同区域）：\n');
fprintf('  真实3麦克风 mean_pd = %.4f\n', mean(pd_real(:)));
fprintf('  贪婪%d节点   mean_pd = %.4f\n', result.n_nodes, result.mean_pd_hist(end));

fprintf('\n演示完成。\n');
end
