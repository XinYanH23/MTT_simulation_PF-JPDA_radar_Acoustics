function ac_cfg = acoustic_config_from_scenario(acousticPos_m, src_height_m)
%ACOUSTIC_CONFIG_FROM_SCENARIO  从场景声学传感器坐标生成 acoustic_config struct。
%
%   ac_cfg = ACOUSTIC_CONFIG_FROM_SCENARIO(acousticPos_m)
%   ac_cfg = ACOUSTIC_CONFIG_FROM_SCENARIO(acousticPos_m, src_height_m)
%
%   acousticPos_m : Mx3 场景世界坐标 [x_east, y_north, z_height] (m)
%                   （与 run_scenario / Step1_Data.mat 一致）
%   src_height_m  : 声源所在平面高度 (m)，默认 median(z_height)+20
%
%   输出 ac_cfg 与 acoustic_config() 字段兼容，但：
%     1) mic_coords 映射为声学场坐标系 [X_east, Y_height, Z_north] (mm)
%     2) 传播参数与 Step2 对齐：SL=110, 球面扩散 n=20, bg=60 dB
%
%   不修改 matlab/acoustic_config.m 原文件。

if nargin < 2 || isempty(src_height_m)
    src_height_m = median(acousticPos_m(:, 3)) + 30.0;
end

M = size(acousticPos_m, 1);

% --- 从 acoustic_config 继承全部默认参数 ---
ac_cfg = acoustic_config();

% --- 坐标：场景 [x, y_north, z_h] (m) → 声学 [X, Y_h, Z_north] (mm) ---
ac_cfg.mic_coords = [acousticPos_m(:,1), acousticPos_m(:,3), acousticPos_m(:,2)] * 1000.0;
ac_cfg.mic_ids    = arrayfun(@(i) num2str(i), 1:M, 'UniformOutput', false);

% --- 传播参数与 Step2 对齐（球面扩散 + 空气吸收）---
%   Step2: recvLevel = SL - 20*log10(r) - alpha*r, snr = recvLevel - bgNoise
%   场模型: Lp = A - n*log10(r) - alpha*r, snr = Lp - noise_floor
ac_cfg.A              = repmat(110.0, M, 1);   % 声源级 SL (dB)
ac_cfg.n              = repmat( 20.0, M, 1);   % 球面扩散
ac_cfg.alpha_atm      = 0.002;                 % dB/m
ac_cfg.noise_floor_db = repmat( 60.0, M, 1);   % 背景噪声 (dB)
ac_cfg.R_ground       = 0.0;                   % 与 Step2 一致：先用自由场

ac_cfg.placement_A_default  = 110.0;
ac_cfg.placement_n_default  =  20.0;
ac_cfg.placement_nf_default =  60.0;

% 检测门限与 Step2 logistic 中点对齐（SNR≈3.5 dB 时 Pd≈0.5）
ac_cfg.det_eta_db     = 3.5;
ac_cfg.det_kappa_db   = 2.0;
ac_cfg.det_threshold  = 0.30;   % 联合 Pd 超过此值即判定 detected（供 PF 融合）
ac_cfg.snr_ref_db     = 3.5;

% --- 声源高度与网格 ---
ac_cfg.src_height_mm = src_height_m * 1000.0;
ac_cfg.placement_area_mm = [0, 1e6, 0, 1e6];   % 1000m×1000m 场景
% 33 站时 10m 网格可接受；若过慢可改为 20000
ac_cfg.grid_res_mm = 15000.0;

end
