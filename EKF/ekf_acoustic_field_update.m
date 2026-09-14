function [x_upd, P_upd, ell_bar] = ekf_acoustic_field_update(x, P, ac_pack, cfg)
%EKF_ACOUSTIC_FIELD_UPDATE  用与 PF 相同的声学似然对高斯状态做乘性再加权。
%
%   PF：  w ← w · ℓ_r · [γ L(x) + (1-γ)]
%   KF：  雷达 PDA 得到 N(x,P) 后，在 σ 点上乘同一 ℓ_a(x)，再匹配均值/协方差。
%
%   ell_bar = E[ℓ_a]（σ 点均值），供 Λ ← Λ_radar · ell_bar，对应
%             ell_fused = ell_radar .* ell_ac 的一阶分解。

nx = numel(x);
x = x(:);
P = (P + P') / 2;
ell_bar = 1;

use_ac = ~isempty(ac_pack) && isstruct(ac_pack) && ...
         isfield(ac_pack, 'detected') && ac_pack.detected;
if ~use_ac
    x_upd = x;
    P_upd = P;
    return
end

[X, w0] = sigma_points_pos(x, P);
if isempty(X)
    ell_bar = acoustic_field_likelihood(x, ac_pack, cfg);
    ell_bar = ell_bar(1);
    x_upd = x;
    P_upd = P;
    return
end

ell = acoustic_field_likelihood(X, ac_pack, cfg);   % 1 × nσ
ell_bar = max(mean(ell), 1e-300);

w = w0 .* ell;
sw = sum(w);
if ~(isfinite(sw) && sw > 1e-300)
    x_upd = x;
    P_upd = P;
    return
end
w = w / sw;

x_upd = X * w(:);
dx = X - x_upd;
P_upd = dx * diag(w) * dx';
P_upd = (P_upd + P_upd') / 2;

% 数值保护：再加权后协方差过小则回退微小过程噪声
mineig = min(eig(P_upd));
if ~(isfinite(mineig) && mineig > 1e-9)
    P_upd = P_upd + 1e-6 * eye(nx);
    P_upd = (P_upd + P_upd') / 2;
end
end

function [X, w0] = sigma_points_pos(x, P)
% 2n+1 对称 σ 点：均值 ± chol 列（1σ）。等权、非负，便于按 ℓ_a 再加权。
% 不用 √(n+κ) 拉伸，避免点落在场外、跳过近处峰值。
nx = numel(x);
try
    L = chol(P, 'lower');
catch
    P = P + 1e-6 * eye(nx);
    [L, p] = chol(P, 'lower');
    if p ~= 0
        X = []; w0 = [];
        return
    end
end
n_sig = 2 * nx + 1;
X = zeros(nx, n_sig);
X(:, 1) = x;
for i = 1:nx
    X(:, 1 + i)      = x + L(:, i);
    X(:, 1 + nx + i) = x - L(:, i);
end
w0 = ones(1, n_sig) / n_sig;
end
