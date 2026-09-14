function [x_fused, P_fused] = imm_fuse(models, mu)
% IMM_FUSE  将各模型估计融合为单一输出。
%
%   [x_fused, P_fused] = imm_fuse(models, mu)
%
%   INPUTS
%     models : 1×M struct（含更新后的 .x .P）
%     mu     : M×1  更新后模型概率
%
%   公式（Bar-Shalom 2001, §11.6.6, eq. 11.6.6-9 ~ 11.6.6-10）：
%
%     x̂  = Σ_j μ_j x̂_j(+)
%
%     P   = Σ_j μ_j [ P_j(+) + (x̂_j − x̂)(x̂_j − x̂)' ]

M = length(models);
mu = mu(:);

nx      = numel(models(1).x);
x_fused = zeros(nx, 1);

for j = 1 : M
    x_fused = x_fused + mu(j) * models(j).x;
end

P_fused = zeros(nx, nx);
for j = 1 : M
    dx      = models(j).x - x_fused;
    P_fused = P_fused + mu(j) * (models(j).P + dx * dx');
end
P_fused = (P_fused + P_fused') / 2;
end
