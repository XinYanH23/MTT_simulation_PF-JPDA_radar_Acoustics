function [pd_joint, xc, zc] = compute_coverage_map(mic_coords, A_vec, n_vec, nf_vec, cfg)
%COMPUTE_COVERAGE_MAP 给定麦克风集合，在覆盖评估网格上计算联合探测概率 P_D_joint(u)。
%
%   [pd_joint, xc, zc] = COMPUTE_COVERAGE_MAP(mic_coords, A_vec, n_vec, nf_vec, cfg)
%
%   输入：
%     mic_coords : Mx3 麦克风世界坐标 [X,Y,Z] (mm)
%     A_vec      : Mx1 各麦克风参考声压级 (dB)
%     n_vec      : Mx1 各麦克风扩散斜率
%     nf_vec     : Mx1 各麦克风环境底噪 (dB)
%   输出：
%     pd_joint   : nx x nz 联合探测概率场 P_D_joint = 1 - prod_m(1-P_D_m)
%     xc, zc     : 覆盖评估网格中心坐标 (mm)
%
%   使用 cfg.placement_area_mm / placement_area_res_mm 定义网格。
%   物理模型与 predict_lp_grid 完全一致（含地面反射），独立实现以支持任意麦克风集合。

M = size(mic_coords, 1);
A_vec  = A_vec(:);
n_vec  = n_vec(:);
nf_vec = nf_vec(:);

% 覆盖评估网格
[xc, zc] = placement_grid(cfg);
[XX, ZZ] = ndgrid(xc, zc);

% 每格 (1 − P_D_m) 的乘积，初始化为 1
one_minus_prod = ones(size(XX));

for m = 1:M
    Lp = lp_grid_for_placement(mic_coords(m,:), A_vec(m), n_vec(m), XX, ZZ, cfg);
    snr = Lp - nf_vec(m);
    pd_m = acoustic_detection_prob(snr, cfg);
    one_minus_prod = one_minus_prod .* (1 - pd_m);
end

pd_joint = 1 - one_minus_prod;

end

% -------------------------------------------------------------------------
function Lp = lp_grid_for_placement(mic, A, n, XX, ZZ, cfg)
%LP_GRID_FOR_PLACEMENT 向量化 SPL，接受显式 mic 参数（不走 cfg.mic_coords 索引）。
dx = XX - mic(1);
dz = ZZ - mic(3);
dy_d = cfg.src_height_mm - mic(2);
r_d  = sqrt(dx.^2 + dy_d.^2 + dz.^2) / 1000.0;
r_d  = max(r_d, 0.05);
Lp_d = A - n .* log10(r_d) - cfg.alpha_atm .* r_d;

if cfg.R_ground > 0
    dy_i = (2 * cfg.ground_y_mm - cfg.src_height_mm) - mic(2);
    r_i  = sqrt(dx.^2 + dy_i.^2 + dz.^2) / 1000.0;
    r_i  = max(r_i, 0.05);
    Lp_i = A - n .* log10(r_i) - cfg.alpha_atm .* r_i;
    I    = 10.^(Lp_d/10) + cfg.R_ground.^2 .* 10.^(Lp_i/10);
    Lp   = 10 * log10(I);
else
    Lp = Lp_d;
end
end

% -------------------------------------------------------------------------
function [xc, zc] = placement_grid(cfg)
%PLACEMENT_GRID 返回覆盖评估网格中心坐标。
if ~isempty(cfg.placement_area_mm)
    ar = cfg.placement_area_mm;
else
    % 自动推断：麦克风包络 + grid_margin_mm
    pts = cfg.mic_coords;
    ar  = [min(pts(:,1)) - cfg.grid_margin_mm, max(pts(:,1)) + cfg.grid_margin_mm, ...
           min(pts(:,3)) - cfg.grid_margin_mm, max(pts(:,3)) + cfg.grid_margin_mm];
end
res = cfg.placement_area_res_mm;
xc  = ar(1):res:ar(2);
zc  = ar(3):res:ar(4);
end
