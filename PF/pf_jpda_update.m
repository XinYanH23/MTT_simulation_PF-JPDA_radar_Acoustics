function [pf, Lambda, ess] = pf_jpda_update(pf, H, ...
                                              Z_valid, R_valid, ...
                                              beta_valid, beta0, lambda_c, ...
                                              ac_pack, cfg)
%PF_JPDA_UPDATE  JPDA 雷达后验混合 + 声学似然修正。
%
%   主路径（cfg.pf_update_mode = 'jpda_posterior'，默认）：
%     1) 预测权 w_- 代表 p^-(x)
%     2) 对每个门内点迹 d：
%           tilde_w_d(m) = w_-(m) * N(z_d; H x^(m), R_d)
%           C_d = sum_m tilde_w_d(m)
%           w_d(m) = tilde_w_d(m) / C_d
%     3) 雷达后验：
%           w_R(m) = beta0 * w_-(m) + sum_d beta_d * w_d(m)
%     4) 声学相对似然 Λ_a(m) = Λ_k^a(p(x^(m)))
%     5) 融合：w(m) ∝ w_R(m) * Λ_a(m)
%
%   旧路径（cfg.pf_update_mode = 'legacy_likelihood'）：
%     w ∝ w_- * (beta0*lambda_c + sum_d beta_d N(z_d;Hx,R)) * ell_a
%     仅用于复现第三章旧版 γ 消融，不再作为论文主算法。

if ~isfield(cfg, 'pf_update_mode') || isempty(cfg.pf_update_mode)
    mode = 'jpda_posterior';
else
    mode = lower(cfg.pf_update_mode);
end

if strcmp(mode, 'legacy_likelihood')
    [pf, Lambda, ess] = pf_jpda_update_legacy(pf, H, ...
        Z_valid, R_valid, beta_valid, beta0, lambda_c, ac_pack, cfg);
    return
end

N     = pf.N;
Dv    = size(Z_valid, 2);
m_dim = size(H, 1);

particles = pf.particles;                 % nx x N
w_minus   = pf.weights(:)';               % 1  x N
w_minus   = max(w_minus, 0);
sw = sum(w_minus);
if sw < 1e-300
    w_minus = ones(1, N) / N;
else
    w_minus = w_minus / sw;
end

%% ── 各关联假设的条件后验 ─────────────────────────────────────────────
C_d   = zeros(1, max(Dv, 1));
w_R   = beta0 * w_minus;

if Dv > 0
    for d = 1:Dv
        z_d   = Z_valid(:, d);
        R_d   = R_valid(:, :, d);
        innov = repmat(z_d, 1, N) - H * particles;

        R_inv    = inv(R_d);
        det_2piR = max((2*pi)^m_dim * det(R_d), 1e-300);
        quad     = sum(innov .* (R_inv * innov), 1);
        gauss    = exp(-0.5 * quad) / sqrt(det_2piR);

        tilde_w = w_minus .* gauss;
        C_d(d)  = max(sum(tilde_w), 1e-300);
        w_d     = tilde_w / C_d(d);
        w_R     = w_R + beta_valid(d) * w_d;
    end
end

w_R = max(w_R, 0);
sum_R = sum(w_R);
if sum_R < 1e-300
    w_R = w_minus;
else
    w_R = w_R / sum_R;
end

%% ── 声学相对似然修正 ─────────────────────────────────────────────────
ell_ac = acoustic_field_likelihood(particles, ac_pack, cfg);

unnorm = w_R .* ell_ac;
sum_w  = sum(unnorm);
if sum_w < 1e-300
    pf.weights = ones(1, N) / N;
else
    pf.weights = unnorm / sum_w;
end

%% ── 模型证据：雷达预测边缘似然 × 声学在雷达后验上的期望 ─────────────
evid_r = beta0;
if Dv > 0
    evid_r = evid_r + sum(beta_valid(:)' .* C_d(1:Dv));
end
evid_a = sum(w_R .* ell_ac);
if ~isfinite(evid_a) || evid_a < 0
    evid_a = 0;
end
Lambda = max(evid_r * evid_a, cfg.pf_lambda_floor);

%% ── ESS + 重采样 ──────────────────────────────────────────────────────
ess = 1 / sum(pf.weights .^ 2);
if ess < cfg.pf_resample_thresh * N
    pf = pf_systematic_resample(pf);
end

pf = pf_extract_state(pf);
end

function [pf, Lambda, ess] = pf_jpda_update_legacy(pf, H, ...
        Z_valid, R_valid, beta_valid, beta0, lambda_c, ac_pack, cfg)
% 旧版 β 加权混合似然，仅用于复现历史实验。

N     = pf.N;
Dv    = size(Z_valid, 2);
m_dim = size(H, 1);

particles = pf.particles;
w_prior   = pf.weights(:)';

ell_radar = beta0 * lambda_c * ones(1, N);
if Dv > 0
    for d = 1:Dv
        z_d   = Z_valid(:, d);
        R_d   = R_valid(:, :, d);
        innov = repmat(z_d, 1, N) - H * particles;
        R_inv    = inv(R_d);
        det_2piR = max((2*pi)^m_dim * det(R_d), 1e-300);
        quad     = sum(innov .* (R_inv * innov), 1);
        gauss    = exp(-0.5 * quad) / sqrt(det_2piR);
        ell_radar = ell_radar + beta_valid(d) * gauss;
    end
end

ell_ac    = acoustic_field_likelihood(particles, ac_pack, cfg);
ell_fused = ell_radar .* ell_ac;
unnorm_w  = w_prior .* ell_fused;
sum_w     = sum(unnorm_w);
if sum_w < 1e-300
    pf.weights = ones(1, N) / N;
else
    pf.weights = unnorm_w / sum_w;
end

Lambda = pf_compute_lambda(ell_fused, cfg.pf_lambda_floor);
ess = 1 / sum(pf.weights .^ 2);
if ess < cfg.pf_resample_thresh * N
    pf = pf_systematic_resample(pf);
end
pf = pf_extract_state(pf);
end
