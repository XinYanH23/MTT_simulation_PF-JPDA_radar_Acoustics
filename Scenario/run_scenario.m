%% =========================================================================
%  Tracker/Scenario/run_scenario.m
%  场景生成入口 —— 城市 / 机场环境 + 无人机异步轨迹
%
%  产出：
%    ① Step1_Data.mat   与 main_imm_jpda_ekf.m 直接对接
%    ② Fig 1  3D 全景（建筑 + 传感器 + 轨迹）
%    ③ Fig 2  俯视图（建筑轮廓 + 路径）
%    ④ Fig 3  各架无人机高度-时间曲线（体现异步出现/消失）
%    ⑤ Fig 4  各架无人机水平速度-时间曲线
%
%  无任何函数重复——所有计算函数仍保留在 v2.0 目录。
%% =========================================================================
clear; clc; close all;
rng(42);

%% ── 路径设置 ────────────────────────────────────────────────────────────
this_dir = fileparts(mfilename('fullpath'));
v2_dir   = fullfile(this_dir, '..', '..', 'v2.0');

addpath(this_dir);
addpath(v2_dir);           % generateUrbanCity / generateUAVTrajectories_v2 等

%% ── 场景选择 ────────────────────────────────────────────────────────────
ScenarioType = 'Urban';    % 'Urban' 或 'Airport'

fprintf('===== Scenario Generator [%s] =====\n', ScenarioType);

%% ── 场景参数 ────────────────────────────────────────────────────────────
switch ScenarioType
    case 'Urban'
        Map.size      = 1000;
        Map.blockSize = 50;
        Map.roadWidth = 20;
        Map.heightMin = 10;
        Map.heightMax = 60;
        UAV.speedRange   = [6, 14];
        UAV.heightMargin = 8;
        AcousticModel.radius   = 150;
        AcousticModel.minRange = 2;
        RadarModel.radius      = 500;
        RadarModel.minRange    = 20;
        buildings = generateUrbanCity(Map);

    case 'Airport'
        Map.size           = 1500;
        Map.centerSize     = 200;
        Map.terminalHeight = 40;
        Map.poleSpacing    = 80;
        Map.poleHeight     = 12;
        UAV.speedRange     = [10, 25];
        UAV.heightMargin   = 15;
        AcousticModel.radius   = 80;
        AcousticModel.minRange = 2;
        RadarModel.radius      = 800;
        RadarModel.minRange    = 20;
        buildings = generateAirportCity(Map);
end

%% ── 通用参数 ────────────────────────────────────────────────────────────
UAV.num = 5;
UAV.dt  = 0.1;
UAV.T   = 50;
time    = (0 : UAV.dt : UAV.T - UAV.dt)';

maxBldgH = max(arrayfun(@(b) b.pos(3) + b.dim(3), buildings));
Map.flyable.xlim = [0, Map.size];
Map.flyable.ylim = [0, Map.size];
Map.flyable.zlim = [maxBldgH + UAV.heightMargin, ...
                    maxBldgH + UAV.heightMargin + 40];
Map.speedRange   = UAV.speedRange;
Map.buildings    = buildings;   % 供碰撞检测使用

%% ── 传感器部署 ──────────────────────────────────────────────────────────
fprintf('正在优化传感器布局...\n');
candidatePos = generateSensorCandidates(buildings, 1.0);
[acousticPos, ~, covScore] = greedySensorPlacement( ...
    candidatePos, Map, AcousticModel, [], 0.99);
radarPos = placeRadarForCoverage(candidatePos, Map, RadarModel, 0.9);

%% ── 无人机轨迹生成 ──────────────────────────────────────────────────────
fprintf('正在生成 %d 架 UAV 异步轨迹...\n', UAV.num);
truth = generateUAVTrajectories_v2(Map, UAV.num);

fprintf('  声学传感器: %d 个  覆盖率: %.1f%%\n', ...
    size(acousticPos,1), covScore*100);
fprintf('  雷达传感器: %d 个\n', size(radarPos,1));
for i = 1:UAV.num
    nv = sum(~isnan(truth(i).pos3D(:,1)));
    fprintf('  UAV #%d: 出生 %.1fs  消失 %.1fs  有效帧 %d\n', ...
        i, truth(i).birthTime, truth(i).deathTime, nv);
end

%% ── 保存（与 main_imm_jpda_ekf.m 接口一致）────────────────────────────
out_mat = fullfile(this_dir, 'Step1_Data.mat');
save(out_mat, 'buildings','acousticPos','radarPos','truth','time','UAV','Map');
fprintf('数据已保存: %s\n\n', out_mat);

%% =========================================================================
%% 可视化
%% =========================================================================
colors_uav = lines(UAV.num);
K          = length(time);

%% ── Figure 1：3D 全景 ───────────────────────────────────────────────────
fig1 = figure('Name','Fig1: 3D Scene Overview', ...
              'Color','w','Position',[50 50 1100 750]);
hold on; grid on; axis equal; box on;
xlabel('X (m)','FontSize',11,'Interpreter','latex');
ylabel('Y (m)','FontSize',11,'Interpreter','latex');
zlabel('Z (m)','FontSize',11,'Interpreter','latex');
title(sprintf('\\bf 3D Scene Overview - %s Scenario', ScenarioType), ...
    'Interpreter','latex','FontSize',14);

% 建筑
for i = 1:length(buildings)
    draw_building_3d(buildings(i));
end

% 声学传感器
if ~isempty(acousticPos)
    scatter3(acousticPos(:,1), acousticPos(:,2), acousticPos(:,3)+1, ...
        60,'b','filled','DisplayName','Acoustic','MarkerEdgeColor','none');
end

% 雷达传感器
if ~isempty(radarPos)
    scatter3(radarPos(:,1), radarPos(:,2), radarPos(:,3)+3, ...
        180,'r','^','filled','DisplayName','Radar','MarkerEdgeColor','k');
end

% UAV 轨迹（彩色线段）
h_legend = gobjects(UAV.num, 1);
for i = 1:UAV.num
    p     = truth(i).pos3D;
    valid = ~isnan(p(:,1));
    if sum(valid) > 2
        h_legend(i) = plot3(p(valid,1), p(valid,2), p(valid,3), ...
            'LineWidth', 2.0, 'Color', colors_uav(i,:), ...
            'DisplayName', sprintf('UAV\\,\\#%d', i));
        % 出生点
        first = find(valid, 1, 'first');
        scatter3(p(first,1), p(first,2), p(first,3), ...
            80, colors_uav(i,:), 'o', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8,'HandleVisibility','off');
        % 消失点
        last = find(valid, 1, 'last');
        scatter3(p(last,1), p(last,2), p(last,3), ...
            80, colors_uav(i,:), 's', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8,'HandleVisibility','off');
    end
end
legend([h_legend; ...
    findobj(fig1,'DisplayName','Acoustic'); ...
    findobj(fig1,'DisplayName','Radar')], ...
    'Location','northeast','Interpreter','latex','FontSize',9);
xlim([-0.05*Map.size, 1.05*Map.size]);
ylim([-0.05*Map.size, 1.05*Map.size]);
zlim([0, Map.flyable.zlim(2)+10]);
view(40, 28);
saveas(fig1, fullfile(this_dir, 'fig1_3d_overview.png'));

%% ── Figure 2：俯视图（XY 平面）─────────────────────────────────────────
fig2 = figure('Name','Fig2: Top-Down View', ...
              'Color','w','Position',[100 100 900 850]);
hold on; grid on; axis equal; box on;
xlabel('X (m)','FontSize',11,'Interpreter','latex');
ylabel('Y (m)','FontSize',11,'Interpreter','latex');
title(sprintf('\\bf Top-Down View - %s Scenario', ScenarioType), ...
    'Interpreter','latex','FontSize',14);

% 建筑轮廓（填充色）
for i = 1:length(buildings)
    b  = buildings(i);
    bx = b.pos(1);  by = b.pos(2);
    bw = b.dim(1);  bd = b.dim(2);
    bh = b.dim(3);
    alpha_val = 0.3 + 0.5 * (bh / (maxBldgH + 1));
    fill([bx, bx+bw, bx+bw, bx], [by, by, by+bd, by+bd], ...
        [0.55 0.55 0.6], 'FaceAlpha', alpha_val, ...
        'EdgeColor',[0.3 0.3 0.35],'LineWidth',0.4,'HandleVisibility','off');
end

% 声学传感器
if ~isempty(acousticPos)
    scatter(acousticPos(:,1), acousticPos(:,2), 55, 'b', 'filled', ...
        'MarkerEdgeColor','w','LineWidth',0.8,'DisplayName','Acoustic');
end

% 雷达传感器（带覆盖圆）
if ~isempty(radarPos)
    for i = 1:size(radarPos,1)
        theta_c = linspace(0, 2*pi, 120);
        xc = radarPos(i,1) + RadarModel.radius * cos(theta_c);
        yc = radarPos(i,2) + RadarModel.radius * sin(theta_c);
        fill(xc, yc, [1 0 0], 'FaceAlpha',0.04,'EdgeColor','r', ...
            'LineWidth',0.8,'HandleVisibility','off');
    end
    scatter(radarPos(:,1), radarPos(:,2), 130, 'r', '^', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.8,'DisplayName','Radar');
end

% 轨迹（用时间渐变色表示进程）
for i = 1:UAV.num
    p     = truth(i).pos3D;
    valid = find(~isnan(p(:,1)));
    if length(valid) < 2, continue; end

    % 渐变：从浅到深
    n_seg = length(valid) - 1;
    base  = colors_uav(i,:);
    for s = 1:n_seg
        alpha_line = 0.3 + 0.7 * s / n_seg;
        col = base * alpha_line + [1 1 1] * (1 - alpha_line);
        col = max(0, min(1, col));
        plot(p(valid(s:s+1), 1), p(valid(s:s+1), 2), ...
            'Color', col, 'LineWidth', 1.6, 'HandleVisibility','off');
    end
    % 图例代理线
    plot(nan, nan, 'Color', base, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('UAV\\,\\#%d', i));
    % 出生 / 消失标记
    scatter(p(valid(1),1), p(valid(1),2), 70, base, 'o', 'filled', ...
        'MarkerEdgeColor','k','HandleVisibility','off');
    scatter(p(valid(end),1), p(valid(end),2), 70, base, 's', 'filled', ...
        'MarkerEdgeColor','k','HandleVisibility','off');
    % 编号标注
    text(p(valid(1),1)+8, p(valid(1),2)+8, sprintf('#%d',i), ...
        'Color',base,'FontSize',9,'FontWeight','bold');
end

legend('Location','northeast','Interpreter','latex','FontSize',9);
xlim([-0.05*Map.size, 1.05*Map.size]);
ylim([-0.05*Map.size, 1.05*Map.size]);
saveas(fig2, fullfile(this_dir, 'fig2_topdown.png'));

%% ── Figure 3：各 UAV 高度-时间曲线（体现异步出现与消失）──────────────
fig3 = figure('Name','Fig3: Altitude Profiles', ...
              'Color','w','Position',[150 150 1000 500]);
hold on; grid on; box on;
xlabel('Time (s)','FontSize',11,'Interpreter','latex');
ylabel('Altitude (m)','FontSize',11,'Interpreter','latex');
title('\bf UAV Altitude Profiles (Asynchronous Birth / Death)', ...
    'Interpreter','latex','FontSize',13);

% 可飞行区间背景
yline(Map.flyable.zlim(1), '--', 'Color',[0.6 0.6 0.6], ...
    'LineWidth',1,'DisplayName','Flight Zone');
yline(Map.flyable.zlim(2), '--', 'Color',[0.6 0.6 0.6], ...
    'LineWidth',1,'HandleVisibility','off');
patch([0 UAV.T UAV.T 0], ...
      [Map.flyable.zlim(1) Map.flyable.zlim(1) ...
       Map.flyable.zlim(2) Map.flyable.zlim(2)], ...
      [0.9 0.95 0.9], 'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');

for i = 1:UAV.num
    p     = truth(i).pos3D;
    valid = ~isnan(p(:,1));
    t_active = time(valid);
    z_active = p(valid, 3);

    plot(t_active, z_active, 'LineWidth', 1.8, 'Color', colors_uav(i,:), ...
        'DisplayName', sprintf('UAV\\,\\#%d (birth %.1fs)', i, truth(i).birthTime));

    % 出生 / 消失竖线
    xline(truth(i).birthTime, ':', 'Color', colors_uav(i,:), ...
        'LineWidth', 1, 'HandleVisibility','off');
    xline(truth(i).deathTime, '-.', 'Color', colors_uav(i,:), ...
        'LineWidth', 1, 'HandleVisibility','off');
end

xlim([0, UAV.T]);
ylim([0, Map.flyable.zlim(2) + 5]);
legend('Location','northeast','Interpreter','latex','FontSize',9);
saveas(fig3, fullfile(this_dir, 'fig3_altitude.png'));

%% ── Figure 4：各 UAV 水平速度-时间曲线 ─────────────────────────────────
fig4 = figure('Name','Fig4: Speed Profiles', ...
              'Color','w','Position',[200 200 1000 500]);
hold on; grid on; box on;
xlabel('Time (s)','FontSize',11,'Interpreter','latex');
ylabel('Horizontal Speed (m/s)','FontSize',11,'Interpreter','latex');
title('\bf UAV Horizontal Speed Profiles', ...
    'Interpreter','latex','FontSize',13);

% 速度范围参考带
yline(UAV.speedRange(1), '--', 'Color',[0.6 0.6 0.6], ...
    'LineWidth',1,'DisplayName','Speed Range');
yline(UAV.speedRange(2), '--', 'Color',[0.6 0.6 0.6], ...
    'LineWidth',1,'HandleVisibility','off');

for i = 1:UAV.num
    if isfield(truth(i),'vel3D')
        v      = truth(i).vel3D;
        valid  = ~isnan(v(:,1));
        t_act  = time(valid);
        spd    = sqrt(v(valid,1).^2 + v(valid,2).^2);
    else
        % 回退：从 pos3D 差分估计速度
        p      = truth(i).pos3D;
        valid  = ~isnan(p(:,1));
        t_act  = time(valid);
        dp     = diff(p(valid,1:2), 1, 1);
        spd    = [0; sqrt(sum(dp.^2, 2))] / UAV.dt;
    end

    plot(t_act, spd, 'LineWidth', 1.8, 'Color', colors_uav(i,:), ...
        'DisplayName', sprintf('UAV\\,\\#%d', i));
end

xlim([0, UAV.T]);
ylim([0, UAV.speedRange(2) * 1.3]);
legend('Location','northeast','Interpreter','latex','FontSize',9);
saveas(fig4, fullfile(this_dir, 'fig4_speed.png'));

fprintf('\n===== 场景生成完成 =====\n');
fprintf('  输出目录: %s\n', this_dir);
fprintf('  图像文件: fig1_3d_overview.png / fig2_topdown.png / fig3_altitude.png / fig4_speed.png\n');
fprintf('  数据文件: Step1_Data.mat\n');

%% =========================================================================
%% 局部辅助函数
%% =========================================================================

function draw_building_3d(b)
% DRAW_BUILDING_3D  在 3D 图中绘制单个建筑（兼容 .pos/.dim 结构）。
if isfield(b,'pos')
    cx = b.pos(1) + b.dim(1)/2;
    cy = b.pos(2) + b.dim(2)/2;
else
    cx = b.center(1);
    cy = b.center(2);
end
w = b.dim(1);   d = b.dim(2);   h = b.dim(3);

% 顶面（灰色）
patch([cx-w/2 cx+w/2 cx+w/2 cx-w/2], ...
      [cy-d/2 cy-d/2 cy+d/2 cy+d/2], ...
      [h h h h], [0.78 0.78 0.82], ...
      'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.4,'FaceAlpha',0.95, ...
      'HandleVisibility','off');

% 侧面（稍暗）
side_col = [0.65 0.65 0.70];
ea       = 0.90;
% 前面（y-）
patch([cx-w/2 cx+w/2 cx+w/2 cx-w/2], ...
      [cy-d/2 cy-d/2 cy-d/2 cy-d/2], [0 0 h h], side_col, ...
      'EdgeColor','none','FaceAlpha',ea,'HandleVisibility','off');
% 右面（x+）
patch([cx+w/2 cx+w/2 cx+w/2 cx+w/2], ...
      [cy-d/2 cy+d/2 cy+d/2 cy-d/2], [0 0 h h], side_col*0.9, ...
      'EdgeColor','none','FaceAlpha',ea,'HandleVisibility','off');
end
