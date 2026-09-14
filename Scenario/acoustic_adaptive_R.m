function R2d = acoustic_adaptive_R(uavPos, sPos, snr_dB, params)
% ACOUSTIC_ADAPTIVE_R  距离/SNR自适应各向异性声学量测噪声协方差（x-y平面）
%
%   物理模型：
%     声学定位误差分为两个方向：
%
%     切向误差（DOA精度）：
%       声学阵列的方位估计标准差满足 Cramér-Rao 界：
%         σ_θ(SNR) = σ_θ0 / sqrt(SNR_linear)
%       弧长误差：σ_t = r · σ_θ(SNR)
%       → 切向误差随距离线性增大，随SNR增大而减小
%
%     径向误差（幅度测距）：
%       声压衰减量给出的距离估计误差：
%         σ_r = σ_r0 · r / sqrt(SNR_linear)
%       → 径向误差通常远大于切向误差（1个量级）
%
%   返回的 R2d 是以传感器-目标连线为主轴的旋转椭圆协方差矩阵。
%
%   INPUTS
%     uavPos  : 1×3 目标位置 [x, y, z]
%     sPos    : 1×3 传感器位置 [x, y, z]
%     snr_dB  : 接收SNR (dB)
%     params  : 参数结构体（可选）
%       .sigma_theta0  基准方位角误差 (rad)，默认 deg2rad(8)
%       .sigma_r0      基准径向误差系数（无量纲），默认 0.15
%       .sigma_min     最小定位标准差 (m)，默认 2.0
%       .sigma_max     最大定位标准差 (m)，默认 50.0
%
%   OUTPUT
%     R2d : 2×2 量测噪声协方差矩阵（世界坐标系 x-y，各向异性椭圆）

if nargin < 4, params = struct(); end

sigma_theta0 = get_param(params, 'sigma_theta0', deg2rad(8));
sigma_r0     = get_param(params, 'sigma_r0',     0.15);
sigma_min    = get_param(params, 'sigma_min',    2.0);
sigma_max    = get_param(params, 'sigma_max',    50.0);

% 传感器到目标的三维距离
dx = uavPos(1) - sPos(1);
dy = uavPos(2) - sPos(2);
r  = sqrt(dx^2 + dy^2 + (uavPos(3) - sPos(3))^2);
r  = max(r, 1.0);

% SNR线性值（下限避免数值奇异）
snr_lin = max(10^(snr_dB / 10), 0.1);

% 各向异性标准差
sigma_t = r * sigma_theta0 / sqrt(snr_lin);   % 切向（弧长）
sigma_r = sigma_r0 * r / sqrt(snr_lin);        % 径向（幅度测距）

% 限幅
sigma_t = min(max(sigma_t, sigma_min), sigma_max);
sigma_r = min(max(sigma_r, sigma_min), sigma_max);

% 旋转矩阵：传感器-目标方位角 → 世界坐标系
phi = atan2(dy, dx);
c   = cos(phi);  s = sin(phi);
Rot = [c, -s; s, c];   % 列向量：[径向方向, 切向方向]

% 局部坐标系协方差（对角，径向×切向）
Sigma_local = diag([sigma_r^2, sigma_t^2]);

% 旋转到世界坐标系
R2d = Rot * Sigma_local * Rot';

end

% ── 内部工具 ──────────────────────────────────────────────────────────────
function v = get_param(s, f, default)
if isfield(s, f)
    v = s.(f);
else
    v = default;
end
end
