function [pf, Lambda, ess, dbg] = pf_jpda_update_debug( ...
    pf, H, Z_valid, R_valid, beta_valid, beta0, lambda_c, ac_pack, cfg)
%PF_JPDA_UPDATE_DEBUG  逐字复现 pf_jpda_update.m，额外返回诊断结构体 dbg。
%  不修改原函数。仅用于 diagnostic_pf_weight_stats.m 诊断脚本。
%
%  额外输出 dbg 字段（均为跨 N 粒子的标量统计量）：
%    w_prior_max/min     更新前权重极值
%    w_unnorm_max/min    乘似然后未归一化权重极值
%    w_norm_max/min      归一化后权重极值
%    log_ell_max/min     log(ell_fused) 极值（对数似然动态范围）
%    ell_radar_max/min/mean  雷达似然极值/均值
%    ell_ac_max/min/mean     声学似然极值/均值（无声学时恒=1）
%    ess_before          更新前 ESS = 1/sum(w_prior²)
%    ess_after           归一化后 ESS = 1/sum(w_norm²)（重采样触发前）
%    resampled           本次是否触发重采样 (0/1)

N     = pf.N;
Dv    = size(Z_valid, 2);
m_dim = size(H, 1);

particles = pf.particles;           % nx x N
w_prior   = pf.weights(:)';         % 1  x N

%% ── 更新前 ESS ────────────────────────────────────────────────────────────
ess_before = 1 / sum(w_prior .^ 2);

%% ── 雷达 PDA 似然（与 pf_jpda_update.m 逐字一致）────────────────────────
ell_radar = beta0 * lambda_c * ones(1, N);

if Dv > 0
    for d = 1:Dv
        z_d   = Z_valid(:, d);
        R_d   = R_valid(:, :, d);
        b_d   = beta_valid(d);
        innov    = repmat(z_d, 1, N) - H * particles;
        R_inv    = inv(R_d);
        det_2piR = (2*pi)^m_dim * det(R_d);
        det_2piR = max(det_2piR, 1e-300);
        quad     = sum(innov .* (R_inv * innov), 1);
        gauss    = exp(-0.5 * quad) / sqrt(det_2piR);
        ell_radar = ell_radar + b_d * gauss;
    end
end

%% ── 声学概率场似然（与 pf_jpda_update / KF 共用）───────────────────────
ell_ac = acoustic_field_likelihood(particles, ac_pack, cfg);

%% ── 融合似然 & 权重更新 ───────────────────────────────────────────────────
ell_fused = ell_radar .* ell_ac;              % 1 x N
unnorm_w  = w_prior   .* ell_fused;           % 1 x N

sum_w = sum(unnorm_w);
if sum_w < 1e-300
    w_norm = ones(1, N) / N;
else
    w_norm = unnorm_w / sum_w;
end
pf.weights = w_norm;

%% ── 模型似然 Λ ────────────────────────────────────────────────────────────
Lambda = pf_compute_lambda(ell_fused, cfg.pf_lambda_floor);

%% ── ESS after（重采样触发前）────────────────────────────────────────────
ess_after = 1 / sum(w_norm .^ 2);
ess = ess_after;

%% ── 重采样（与原函数一致）────────────────────────────────────────────────
did_resample = 0;
if ess_after < cfg.pf_resample_thresh * N
    pf = pf_systematic_resample(pf);
    did_resample = 1;
end

pf = pf_extract_state(pf);

%% ── 构造 dbg 结构体 ───────────────────────────────────────────────────────
log_ell = log(max(ell_fused, 1e-300));   % 避免 log(0)

dbg.w_prior_max    = max(w_prior);
dbg.w_prior_min    = min(w_prior);
dbg.w_unnorm_max   = max(unnorm_w);
dbg.w_unnorm_min   = min(unnorm_w);
dbg.w_norm_max     = max(w_norm);
dbg.w_norm_min     = min(w_norm);
dbg.log_ell_max    = max(log_ell);
dbg.log_ell_min    = min(log_ell);
dbg.ell_radar_max  = max(ell_radar);
dbg.ell_radar_min  = min(ell_radar);
dbg.ell_radar_mean = mean(ell_radar);
dbg.ell_ac_max     = max(ell_ac);
dbg.ell_ac_min     = min(ell_ac);
dbg.ell_ac_mean    = mean(ell_ac);
dbg.ess_before     = ess_before;
dbg.ess_after      = ess_after;
dbg.resampled      = did_resample;
end
