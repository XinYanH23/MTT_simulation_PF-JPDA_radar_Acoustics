%% =========================================================
%  Step2 v2.5 : 异构物理探测仿真（雷达 + 声学）
%
%  修改记录（v2.5）：
%    - 修复：目标速度改用 truth.vel3D 真值速度
%      （原代码硬编码 [10 5 0]，所有目标速度相同）
%    - 修复：vel3D 缺失时自动用有限差分估计速度（向后兼容旧 mat）
%    - 维持：声学物理模型、遮挡检测、objectDetection 封装不变
%% =========================================================

fprintf('===== Step2 v2.5: 异构探测仿真启动 =====\n');

%% 0. 加载 Step1（若工作区无场景变量）
if ~exist('truth', 'var') || ~exist('radarPos', 'var')
    s1 = fullfile(fileparts(mfilename('fullpath')), 'Step1_Data.mat');
    if ~exist(s1, 'file')
        error('缺少 Step1_Data.mat，请先运行 run_scenario.m');
    end
    load(s1);
    fprintf('  已加载: %s\n', s1);
end

%% 1. 参数配置
Env.bgNoise         = 60;    % dB，背景噪声 SPL
AcParam.SL          = 110;   % dB，无人机声源级
AcParam.alpha       = 0.002; % dB/m，空气吸收系数
AcParam.threshold   = 3.5;   % dB，检测门限（SNR）

% 自适应各向异性噪声参数（替代原固定 sigmaPos）
AcParam.R.sigma_theta0 = deg2rad(8);  % 基准DOA角度误差（Cramér-Rao界）
AcParam.R.sigma_r0     = 0.15;        % 基准幅度测距误差系数
AcParam.R.sigma_min    = 2.0;         % 最小定位标准差 (m)
AcParam.R.sigma_max    = 50.0;        % 最大定位标准差 (m)
AcParam.sigma_z        = 10.0;        % 垂直方向固定误差 (m)，声学高度估计差
RadarParam.ReferenceRCS = -12;  % dBsm，参考 RCS

%% 2. 初始化雷达探测器（Radar Toolbox）
numRadars = size(radarPos, 1);
radars    = cell(numRadars, 1);
for i = 1 : numRadars
    radars{i} = radarDataGenerator( ...
        'SensorIndex',       100 + i,       ...
        'UpdateRate',        1 / UAV.dt,    ...
        'MountingLocation',  radarPos(i,:), ...
        'FieldOfView',       [120, 60],     ...
        'RangeLimits',       [20, 1000],    ...
        'RangeResolution',   1.0,           ...
        'AzimuthResolution', 1.0,           ...
        'ReferenceRCS',      RadarParam.ReferenceRCS, ...
        'HasElevation',      true,          ...
        'FalseAlarmRate',    1e-6);
end

%% 3. 主仿真循环
numSteps       = length(time);
allDetections  = cell(numSteps, 1);

for t = 1 : numSteps
    currTime  = time(t);
    frameDets = [];

    %% ── A. 构造目标平台信息 ──────────────────────────────────────────
    targetPlatforms = struct( ...
        'PlatformID', {}, 'Position', {}, 'Velocity', {}, ...
        'Orientation', {}, 'Size', {});
    count = 0;

    for n = 1 : length(truth)
        % 检查当前帧是否有有效位置
        if t > size(truth(n).pos3D, 1), continue; end
        pos_t = truth(n).pos3D(t, :);
        if any(isnan(pos_t)), continue; end

        count = count + 1;
        targetPlatforms(count).PlatformID  = n;
        targetPlatforms(count).Position    = pos_t;

        %% ── 速度：优先使用真值速度 vel3D ──────────────────────────
        if isfield(truth(n), 'vel3D') && ...
                t <= size(truth(n).vel3D, 1) && ...
                ~any(isnan(truth(n).vel3D(t, :)))
            % v2.2+ 的 truth 含有 vel3D 字段
            vel_t = truth(n).vel3D(t, :);
        else
            % 向后兼容：旧 truth 无 vel3D，用有限差分估计
            if t < size(truth(n).pos3D, 1) && ~any(isnan(truth(n).pos3D(t+1,:)))
                vel_t = (truth(n).pos3D(t+1,:) - pos_t) / UAV.dt;
            elseif t > 1 && ~any(isnan(truth(n).pos3D(t-1,:)))
                vel_t = (pos_t - truth(n).pos3D(t-1,:)) / UAV.dt;
            else
                vel_t = [0, 0, 0];   % 边界帧：速度置零
            end
        end

        targetPlatforms(count).Velocity    = vel_t;
        targetPlatforms(count).Orientation = quaternion([1, 0, 0, 0]);
        targetPlatforms(count).Size        = struct('Length',0.5,'Width',0.5,'Height',0.2);
    end

    if isempty(targetPlatforms), continue; end

    %% ── B. 雷达探测 ─────────────────────────────────────────────────
    for r = 1 : numRadars
        [rdDets, ~, ~] = radars{r}(targetPlatforms, currTime);
        if isempty(rdDets), continue; end

        for k = 1 : length(rdDets)
            meas3d = rdDets{k}.Measurement';     % 转为行向量以供遮挡检测
            if ~isOccluded(meas3d, radarPos(r,:), buildings)
                rdDets{k}.ObjectAttributes = struct( ...
                    'ClassConfidence', 0.65, ...
                    'SensorType',      'Radar');
                frameDets = [frameDets; rdDets{k}];  %#ok<AGROW>
            end
        end
    end

    %% ── C. 声学探测 ─────────────────────────────────────────────────
    for n = 1 : length(targetPlatforms)
        uavPos = targetPlatforms(n).Position;   % 1×3

        for s = 1 : size(acousticPos, 1)
            sPos = acousticPos(s, :);
            dist = norm(uavPos - sPos);

            % 球面扩散 + 空气吸收衰减（dB 域）
            recvLevel = AcParam.SL ...
                      - 20 * log10(max(dist, 1)) ...
                      - AcParam.alpha * dist;
            snr = recvLevel - Env.bgNoise;

            % Logistic 检测概率模型（随SNR连续变化）
            Pd = 1 / (1 + exp(-(snr - AcParam.threshold)));

            if rand < Pd
                if ~isOccluded(uavPos, sPos, buildings)
                    % ── 自适应各向异性量测噪声 ────────────────────────
                    % R2d 为以传感器-目标连线为主轴的旋转椭圆：
                    %   径向（幅度测距）误差 >> 切向（DOA角度）误差
                    %   两者均随距离增大、随SNR降低而增大
                    R2d = acoustic_adaptive_R(uavPos, sPos, snr, AcParam.R);

                    % 从各向异性分布采样水平噪声（chol 保证实数采样）
                    R2d_sym = (R2d + R2d') / 2;
                    L       = chol(R2d_sym, 'lower');
                    noise_xy = (L * randn(2, 1))';
                    noise_z  = AcParam.sigma_z * randn;
                    measPos  = uavPos + [noise_xy, noise_z];

                    % 3×3 量测噪声（z固定，xy为自适应椭圆）
                    R3d = blkdiag(R2d, AcParam.sigma_z^2);

                    acDet = objectDetection(currTime, measPos', ...
                        'SensorIndex', s, ...
                        'MeasurementNoise', R3d, ...
                        'ObjectClassID', 1);

                    acDet.ObjectAttributes = struct( ...
                        'ClassConfidence', max(0.3, min(0.95, Pd)), ...
                        'SensorType',      'Acoustic', ...
                        'SNR_dB',          snr);

                    frameDets = [frameDets; acDet];  %#ok<AGROW>
                end
            end
        end
    end

    allDetections{t} = frameDets;

    if mod(t, 100) == 0
        n_det = 0;
        if ~isempty(frameDets), n_det = length(frameDets); end
        fprintf('  帧 %4d/%d | 活跃目标: %d | 当前探测: %d\n', ...
            t, numSteps, length(targetPlatforms), n_det);
    end
end

%% 4. 保存
save('Step2_HeteroDetections.mat', 'allDetections', 'time', 'truth', 'buildings');
fprintf('===== Step2 完成 | 数据已保存至 Step2_HeteroDetections.mat =====\n');

%% ──────────────────────────────────────────────────────────────────────
%% 局部辅助函数（必须放在脚本末尾）

function occluded = isOccluded(p1, p2, buildings)
% ISOCCLUDED  判断 p1→p2 连线是否被任意建筑体遮挡。
%   p1, p2 均为 1×3 行向量。
%   建筑须含 .pos [x y z] 和 .dim [w d h] 字段。
occluded = false;
for i = 1 : length(buildings)
    b   = buildings(i);
    box = [b.pos(1),            b.pos(1) + b.dim(1), ...
           b.pos(2),            b.pos(2) + b.dim(2), ...
           b.pos(3),            b.pos(3) + b.dim(3)];
    if lineBoxIntersection(p1, p2, box)
        occluded = true;
        return
    end
end
end

function hit = lineBoxIntersection(p1, p2, box)
% LINEBOXINTERSECTION  参数化线段与 AABB 包围盒的快速相交测试。
%   box = [xmin xmax ymin ymax zmin zmax]
d    = p2 - p1;
tmin = -inf;
tmax =  inf;
for i = 1 : 3
    if abs(d(i)) > 1e-12
        t1 = (box(2*i-1) - p1(i)) / d(i);
        t2 = (box(2*i)   - p1(i)) / d(i);
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));
    elseif p1(i) < box(2*i-1) || p1(i) > box(2*i)
        hit = false;
        return
    end
end
hit = (tmax >= tmin) && (tmax >= 0) && (tmin <= 1);
end
