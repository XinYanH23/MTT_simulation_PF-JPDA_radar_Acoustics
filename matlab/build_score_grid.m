function [S, xc, zc, XX, ZZ] = build_score_grid(measured_lp, snr_db, w, cfg)
%BUILD_SCORE_GRID 多节点 SPL 协同：在 X-Z 网格上构造融合对数似然得分场 S(u)。
%
%   [S, xc, zc, XX, ZZ] = BUILD_SCORE_GRID(measured_lp, snr_db, w, cfg)
%
%   每个网格点 u 处：
%     S(u) = -0.5 * sum_m  w_m * ( (Lp_meas_m - Lp_pred_m(u)) / sigma_m )^2
%   其中：
%     Lp_pred_m(u) 由 predict_lp_grid 给出（含地面反射）
%     sigma_m = sigma0 * exp(beta * max(0, snr_ref - snr_m))   异方差膨胀
%     w_m      = 节点级 SNR 可靠性权重（node_snr_weights）
%
%   S(u) 是未归一化对数似然得分（越大越匹配），后续由三种映射转为概率场。
%   网格范围由麦克风包络 + margin 决定，分辨率 cfg.grid_res_mm。

% 网格范围（来自麦克风包络）
pts = cfg.mic_coords;
xmin = min(pts(:,1)) - cfg.grid_margin_mm;
xmax = max(pts(:,1)) + cfg.grid_margin_mm;
zmin = min(pts(:,3)) - cfg.grid_margin_mm;
zmax = max(pts(:,3)) + cfg.grid_margin_mm;

xc = xmin:cfg.grid_res_mm:xmax;
zc = zmin:cfg.grid_res_mm:zmax;
[XX, ZZ] = ndgrid(xc, zc);

measured_lp = measured_lp(:);
snr_db = snr_db(:);
w = w(:);
M = numel(measured_lp);

chi2 = zeros(size(XX));
for m = 1:M
    if isnan(measured_lp(m)) || w(m) <= 0
        continue;   % 跳过缺测/零权重节点
    end
    Lp_pred = predict_lp_grid(m, XX, ZZ, cfg.src_height_mm, cfg);
    sigma_m = cfg.sigma0_db * exp(cfg.beta_snr * max(0, cfg.snr_ref_db - snr_db(m)));
    r = (measured_lp(m) - Lp_pred) / sigma_m;
    chi2 = chi2 + w(m) * r.^2;
end

S = -0.5 * chi2;

end
