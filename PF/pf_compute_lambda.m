function Lambda = pf_compute_lambda(ell_unnorm, lambda_floor)
%PF_COMPUTE_LAMBDA  粒子似然均值，供 imm_update_mu 使用。
%
%   ell_unnorm : 1 x N，更新后未归一化的原始似然值 ell(m)
%                （注意：不是乘以 1/N 的 w_prior，是 ell_fused 本身）
%   Lambda     = mean(ell_unnorm)，与各模型之间的相对大小一致

if nargin < 2 || isempty(lambda_floor); lambda_floor = 1e-300; end
Lambda = max(mean(ell_unnorm(:)), lambda_floor);
end
