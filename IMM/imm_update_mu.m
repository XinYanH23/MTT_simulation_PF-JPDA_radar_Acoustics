function mu_new = imm_update_mu(mu, Lambda, c_bar)
% IMM_UPDATE_MU  根据各模型的似然更新模型概率。
%
%   mu_new = imm_update_mu(mu, Lambda, c_bar)
%
%   INPUTS
%     mu     : M×1  当前模型概率（混合前）
%     Lambda : M×1  各模型似然（PDA 归一化常数，见 ekf_pda_update）
%     c_bar  : M×1  imm_mix 返回的预测模型概率
%
%   公式（Bar-Shalom 2001, §11.6.6, eq. 11.6.6-8）：
%
%     μ_j(+) = Λ_j · c̄_j / Σ_k ( Λ_k · c̄_k )
%
%   物理含义：
%     c̄_j 是从模型转移矩阵得到的先验概率；
%     Λ_j 是观测对模型 j 的支持度（似然）；
%     乘积再归一化即贝叶斯后验模型概率。

Lambda = Lambda(:);
c_bar  = c_bar(:);

unnorm = Lambda .* c_bar;
denom  = sum(unnorm);

if denom < 1e-300
    mu_new = ones(length(mu), 1) / length(mu);     % 退化：均匀分布
else
    mu_new = unnorm / denom;
end

% 数值保护
mu_new = max(mu_new, 1e-6);
mu_new = mu_new / sum(mu_new);
end
