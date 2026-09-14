function [Z, R_list, sensor_types] = detections_to_ZR(detections_cell, cfg)
% DETECTIONS_TO_ZR  将 Step2 输出的 allDetections{k} 转换为
%                    标准量测矩阵 Z 和噪声列表 R_list。
%
%   [Z, R_list, sensor_types] = detections_to_ZR(detections_cell, cfg)
%
%   INPUTS
%     detections_cell : allDetections{k}，来自 Step2_HeteroDetections.mat
%                       每个元素为 objectDetection 对象或 struct，含：
%                         .Measurement         — 位置量测向量 (≥2×1)
%                         .MeasurementNoise    — 量测噪声协方差 (≥2×2)
%                         .ObjectAttributes    — 含 SensorType 字段
%     cfg             : config_tracker 输出
%
%   OUTPUTS
%     Z           : 2×D  位置量测矩阵（只取 x, y 分量）
%     R_list      : 2×2×D  对应量测噪声
%     sensor_types: 1×D cell  各探测的传感器类型字符串
%
%   接口设计使跟踪器与传感器类型解耦：
%     雷达探测  → 使用 cfg.R_radar（除非 MeasurementNoise 更精确）
%     声学探测  → 使用 cfg.R_acoustic
%     其他/未知 → 使用 cfg.R_acoustic（保守估计）

if isempty(detections_cell)
    Z            = zeros(2, 0);
    R_list       = zeros(2, 2, 0);
    sensor_types = {};
    return
end

% 处理 cell 数组或 struct 数组
if iscell(detections_cell)
    det_list = detections_cell;
    n_det    = length(det_list);
else
    n_det    = length(detections_cell);
    det_list = num2cell(detections_cell);
end

Z            = zeros(2, n_det);
R_list       = zeros(2, 2, n_det);
sensor_types = cell(1, n_det);
count        = 0;

for i = 1 : n_det
    det = det_list{i};

    %% 提取量测值
    if isobject(det)
        meas  = det.Measurement(:);
        noise = det.MeasurementNoise;
        attr  = det.ObjectAttributes;
    elseif isstruct(det)
        meas  = det.Measurement(:);
        noise = det.MeasurementNoise;
        attr  = det.ObjectAttributes;
    else
        continue
    end

    % 只保留 x, y 分量（忽略高度 z 等）
    if length(meas) < 2, continue; end
    z2d = meas(1:2);

    %% 确定传感器类型与量测噪声
    sensor_type = 'Unknown';
    if iscell(attr)
        for ai = 1:length(attr)
            a = attr{ai};
            if isstruct(a) && isfield(a, 'SensorType')
                sensor_type = a.SensorType;
                break
            end
        end
    elseif isstruct(attr) && isfield(attr, 'SensorType')
        sensor_type = attr.SensorType;
    end

    % 优先使用量测本身携带的噪声；若不可信则用默认值
    if size(noise, 1) >= 2 && size(noise, 2) >= 2
        R2d = noise(1:2, 1:2);
        % 检查正定性
        if ~is_pos_def(R2d)
            R2d = default_R(sensor_type, cfg);
        end
    else
        R2d = default_R(sensor_type, cfg);
    end

    count = count + 1;
    Z(:, count)          = z2d;
    R_list(:, :, count)  = R2d;
    sensor_types{count}  = sensor_type;
end

% 截断为实际探测数
Z            = Z(:, 1:count);
R_list       = R_list(:, :, 1:count);
sensor_types = sensor_types(1:count);
end

%% ── 辅助函数 ──────────────────────────────────────────────────────────
function R = default_R(sensor_type, cfg)
if strcmpi(sensor_type, 'Radar')
    R = cfg.R_radar;
else
    R = cfg.R_acoustic;
end
end

function tf = is_pos_def(M)
try
    chol(M);
    tf = true;
catch
    tf = false;
end
end
