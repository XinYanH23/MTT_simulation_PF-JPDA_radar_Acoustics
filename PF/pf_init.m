function pf = pf_init(x0, P0, N, cov_scale)
%PF_INIT  从高斯初始分布采样 N 粒子，初始化单模型 PF struct。
%
%   pf = PF_INIT(x0, P0, N)
%   pf = PF_INIT(x0, P0, N, cov_scale)
%
%   输出 pf struct 字段：
%     .particles  nx x N   粒子矩阵（每列一个粒子）
%     .weights    1  x N   归一化权重（初始均匀）
%     .N          粒子数
%     .x          nx x 1   加权均值（与 KF model.x 对齐）
%     .P          nx x nx  加权协方差（与 KF model.P 对齐）

if nargin < 4 || isempty(cov_scale); cov_scale = 1.0; end
nx = numel(x0);
x0 = x0(:);
L  = chol(P0 * cov_scale, 'lower');
particles = repmat(x0, 1, N) + L * randn(nx, N);
pf = struct('particles', particles, ...
            'weights',   ones(1,N)/N, ...
            'N',         N, ...
            'x',         x0, ...
            'P',         P0);
end
