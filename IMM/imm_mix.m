function [models_mixed, c_bar] = imm_mix(models, mu, Pi)
% IMM_MIX  IMM 交互/混合步骤。
%
%   [models_mixed, c_bar] = imm_mix(models, mu, Pi)
%
%   计算每个模型 j 的混合初始条件 x̄_{0j}、P̄_{0j}，
%   供后续 EKF 预测步骤使用。
%
%   INPUTS
%     models : 1×M struct（含 .x .P）
%     mu     : M×1  当前模型概率
%     Pi     : M×M  Markov 转移矩阵，Pi(i,j) = P(j|i)
%
%   OUTPUTS
%     models_mixed : 1×M struct（混合后 .x .P 已更新）
%     c_bar        : M×1  预测模型概率 c̄_j = Σ_i π_{ij} μ_i
%
%   公式（Bar-Shalom 2001, §11.6.6, eq. 11.6.6-2 ~ 11.6.6-4）：
%
%     c̄_j       = Σ_i π_{ij} μ_i
%
%     μ_{i|j}   = π_{ij} μ_i / c̄_j        （混合系数）
%
%     x̄_{0j}   = Σ_i μ_{i|j} x̂_i          （混合状态）
%
%     P̄_{0j}   = Σ_i μ_{i|j} [ P_i + Δx_i Δx_i' ]
%               where Δx_i = x̂_i − x̄_{0j}  （混合协方差）

M            = length(models);
mu           = mu(:);
models_mixed = models;              % 复制结构，下面逐模型覆盖

% 预测模型概率
c_bar = Pi' * mu;                   % M×1
c_bar = max(c_bar, 1e-300);         % 防止除零

for j = 1 : M
    % 混合系数 μ_{i|j}
    mu_ij = Pi(:, j) .* mu ./ c_bar(j);   % M×1

    % 混合状态
    nx    = numel(models(1).x);
    x_mix = zeros(nx, 1);
    for i = 1 : M
        x_mix = x_mix + mu_ij(i) * models(i).x;
    end

    % 混合协方差（含均值散度项）
    P_mix = zeros(nx, nx);
    for i = 1 : M
        dx    = models(i).x - x_mix;
        P_mix = P_mix + mu_ij(i) * (models(i).P + dx * dx');
    end
    P_mix = (P_mix + P_mix') / 2;      % 强制对称

    models_mixed(j).x = x_mix;
    models_mixed(j).P = P_mix;
end
end
