function ell_ac = acoustic_field_likelihood(X, ac_pack, cfg)
%ACOUSTIC_FIELD_LIKELIHOOD  在状态位置上查询当前帧声学相对似然。
%
%   主路径：Λ_a(u) = exp[S(u) - max S]，即声学空间相对似然场。
%   本帧声学无效（未检测）时返回 1，算法退化为雷达后验。
%
%   旧路径（cfg.pf_update_mode = 'legacy_likelihood'）：
%     查询 Softmax 场，再做 γ L + (1-γ)。仅用于复现历史实验。

N = size(X, 2);
ell_ac = ones(1, N);

use_ac = ~isempty(ac_pack) && isstruct(ac_pack) && ...
         isfield(ac_pack, 'detected') && ac_pack.detected;
if ~use_ac
    return
end

legacy = isfield(cfg, 'pf_update_mode') && ...
         strcmpi(cfg.pf_update_mode, 'legacy_likelihood');

kappa = 1000;
if isfield(cfg, 'pf_m2ac_scale') && ~isempty(cfg.pf_m2ac_scale)
    kappa = cfg.pf_m2ac_scale;
end
bg = 1e-6;
if isfield(cfg, 'pf_ac_bg_likelihood') && ~isempty(cfg.pf_ac_bg_likelihood)
    bg = cfg.pf_ac_bg_likelihood;
end

px_mm = X(1, :) * kappa;
py_mm = X(2, :) * kappa;

if ~isfield(ac_pack, 'xc') || ~isfield(ac_pack, 'zc') ...
        || isempty(ac_pack.xc) || isempty(ac_pack.zc)
    return
end
xc = ac_pack.xc(:)';
zc = ac_pack.zc(:)';

if legacy
    L_field = legacy_softmax_field(ac_pack);
else
    L_field = relative_likelihood_field(ac_pack);
end
if isempty(L_field)
    return
end

try
    F_interp = griddedInterpolant({xc, zc}, double(L_field), 'linear', 'none');
    ell_interp = F_interp(px_mm(:), py_mm(:));
    ell_interp = ell_interp(:)';
    ell_interp(isnan(ell_interp)) = bg;
catch
    ell_interp = bg * ones(1, N);
end
ell_interp = max(ell_interp, bg);

if legacy
    if isfield(cfg, 'pf_ac_gamma_mode') && strcmpi(cfg.pf_ac_gamma_mode, 'snr')
        gamma = min(max(ac_pack.p_detect_joint, 0), 1);
    else
        gamma = 0.3;
        if isfield(cfg, 'pf_ac_gamma_fixed') && ~isempty(cfg.pf_ac_gamma_fixed)
            gamma = cfg.pf_ac_gamma_fixed;
        end
    end
    ell_ac = gamma * ell_interp + (1 - gamma);
else
    ell_ac = ell_interp;
end
end

function L_field = relative_likelihood_field(ac_pack)
if isfield(ac_pack, 'L_rel') && ~isempty(ac_pack.L_rel)
    L_field = ac_pack.L_rel;
    return
end
if isfield(ac_pack, 'S') && ~isempty(ac_pack.S)
    S = double(ac_pack.S);
    Smax = max(S(:));
    if ~isfinite(Smax)
        L_field = [];
        return
    end
    L_field = exp(S - Smax);
    return
end
L_field = [];
end

function L_field = legacy_softmax_field(ac_pack)
L_field = [];
method_ac = 'softmax';
if isfield(ac_pack, 'method') && ~isempty(ac_pack.method)
    method_ac = lower(ac_pack.method);
end
switch method_ac
    case 'intensity'
        if isfield(ac_pack, 'L_intensity'); L_field = ac_pack.L_intensity; end
    case 'nms_gmm'
        if isfield(ac_pack, 'L_nms_gmm'); L_field = ac_pack.L_nms_gmm; end
    otherwise
        if isfield(ac_pack, 'L_softmax'); L_field = ac_pack.L_softmax; end
end
end
