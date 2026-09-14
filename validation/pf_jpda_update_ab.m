function [pf, Lambda, ess, dbg] = pf_jpda_update_ab(pf, H, ...
    Z_valid, R_valid, beta_valid, beta0, lambda_c, ac_pack, cfg, ac_field_mode)
%PF_JPDA_UPDATE_AB  A/B 验证专用：在 post 与 L_softmax 间切换声学场。
%  不修改 PF/pf_jpda_update.m；仅用于 validation/run_acoustic_likelihood_ab.m。
%
%   ac_field_mode : 'post' | 'L_softmax' | 'L_intensity' | 'L_nms_gmm'

if nargin < 10 || isempty(ac_field_mode)
    ac_field_mode = 'L_softmax';
end

N     = pf.N;
Dv    = size(Z_valid, 2);
m_dim = size(H, 1);

particles = pf.particles;
w_prior   = pf.weights(:)';

ess_before = 1 / sum(w_prior .^ 2);

%% 雷达 PDA 似然（与 pf_jpda_update.m 一致）
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

%% 声学场似然
ell_ac = ones(1, N);
use_ac = ~isempty(ac_pack) && isstruct(ac_pack) && ac_pack.detected;

if use_ac
    kappa  = cfg.pf_m2ac_scale;
    px_mm  = particles(1, :) * kappa;
    py_mm  = particles(2, :) * kappa;
    xc     = ac_pack.xc(:)';
    zc     = ac_pack.zc(:)';

    mode_l = lower(ac_field_mode);
    switch mode_l
        case 'post'
            if ~isfield(ac_pack, 'post') || isempty(ac_pack.post)
                error('pf_jpda_update_ab:post_missing', 'ac_pack.post 不可用');
            end
            L_ac_field = ac_pack.post;
        case 'intensity'
            L_ac_field = ac_pack.L_intensity;
        case 'nms_gmm'
            L_ac_field = ac_pack.L_nms_gmm;
        otherwise
            if ~isfield(ac_pack, 'L_softmax') || isempty(ac_pack.L_softmax)
                error('pf_jpda_update_ab:L_softmax_missing', 'ac_pack.L_softmax 不可用');
            end
            L_ac_field = ac_pack.L_softmax;
    end

    bg = cfg.pf_ac_bg_likelihood;
    try
        F_interp   = griddedInterpolant({xc, zc}, double(L_ac_field), 'linear', 'none');
        ell_interp = F_interp(px_mm(:), py_mm(:));
        ell_interp = ell_interp(:)';
        ell_interp(isnan(ell_interp)) = bg;
    catch
        ell_interp = bg * ones(1, N);
    end
    ell_interp = max(ell_interp, bg);

    if strcmpi(cfg.pf_ac_gamma_mode, 'snr')
        gamma = min(max(ac_pack.p_detect_joint, 0), 1);
    else
        gamma = cfg.pf_ac_gamma_fixed;
    end
    ell_ac = gamma * ell_interp + (1 - gamma);
end

ell_fused = ell_radar .* ell_ac;
unnorm_w  = w_prior .* ell_fused;

sum_w = sum(unnorm_w);
if sum_w < 1e-300
    w_norm = ones(1, N) / N;
else
    w_norm = unnorm_w / sum_w;
end
pf.weights = w_norm;

Lambda = pf_compute_lambda(ell_fused, cfg.pf_lambda_floor);

ess_after = 1 / sum(w_norm .^ 2);
ess = ess_after;

did_resample = 0;
if ess_after < cfg.pf_resample_thresh * N
    pf = pf_systematic_resample(pf);
    did_resample = 1;
end

pf = pf_extract_state(pf);

if nargout >= 4
    log_ell = log(max(ell_fused, 1e-300));
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
    dbg.ac_field_mode  = ac_field_mode;
end
end
