function Lp = predict_lp_grid(mic_idx, XX, ZZ, height_mm, cfg)
%PREDICT_LP_GRID 向量化版单点 SPL 模型：在整个 X-Z 网格上预测某麦克风的声压级。
%
%   Lp = PREDICT_LP_GRID(mic_idx, XX, ZZ, height_mm, cfg)
%
%   输入：
%     mic_idx   : 麦克风行号
%     XX, ZZ    : ndgrid 生成的网格坐标 (mm)，同尺寸 nx x nz
%     height_mm : 声源高度 Y (mm)，标量
%     cfg       : 参数 struct
%   输出：
%     Lp : nx x nz 预测声压级 (dB)
%
%   与 single_mic_spl 完全一致的物理模型，仅做向量化以加速建网格。

mic = cfg.mic_coords(mic_idx, :);
A = cfg.A(mic_idx);
n = cfg.n(mic_idx);

dx = XX - mic(1);
dz = ZZ - mic(3);

% 直达路径
dy_d = height_mm - mic(2);
r_d = sqrt(dx.^2 + dy_d.^2 + dz.^2) / 1000.0;
r_d = max(r_d, 0.05);
Lp_d = A - n .* log10(r_d) - cfg.alpha_atm .* r_d;

if cfg.R_ground > 0
    % 镜像源高度
    dy_i = (2 * cfg.ground_y_mm - height_mm) - mic(2);
    r_i = sqrt(dx.^2 + dy_i.^2 + dz.^2) / 1000.0;
    r_i = max(r_i, 0.05);
    Lp_i = A - n .* log10(r_i) - cfg.alpha_atm .* r_i;
    I = 10.^(Lp_d / 10) + cfg.R_ground.^2 .* 10.^(Lp_i / 10);
    Lp = 10 * log10(I);
else
    Lp = Lp_d;
end

end
