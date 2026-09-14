function [models_mixed, c_bar] = imm_pf_mix(models, mu, Pi)
%IMM_PF_MIX  IMM 粒子混合（替代 KF 版 imm_mix）。
%
%   [models_mixed, c_bar] = IMM_PF_MIX(models, mu, Pi)
%
%   INPUTS
%     models : 1×M struct，每元素含 .particles (nx×N), .weights (1×N), .N
%     mu     : M×1 当前模型概率
%     Pi     : M×M Markov 转移矩阵
%
%   OUTPUTS
%     models_mixed : 混合后粒子集（.particles, .weights 更新，均值/协方差同步）
%     c_bar        : M×1 预测模型概率（与 imm_mix 接口完全一致）
%
%   ── 粒子混合算法 ──────────────────────────────────────────────────
%   对模型 j，混合后分布：p_{0j}(x) = sum_i mu_{i|j} * p_i(x)
%
%   对每个粒子 m = 1..N：
%     1. 以概率 mu_{i|j} 选择源模型 i*
%     2. 从模型 i* 按 w_{i*} 重采样一个粒子
%     => 等价于从混合分布中 Monte Carlo 采样

M  = length(models);
mu = mu(:);
N  = models(1).N;
nx = size(models(1).particles, 1);

c_bar = Pi' * mu;
c_bar = max(c_bar, 1e-300);

% 预先建立各模型的累积权重，避免内层循环重算
CW = cell(M, 1);
for i = 1:M
    w = max(models(i).weights(:), 0);
    w = w / sum(w);
    CW{i} = [0; cumsum(w(:))];   % (N+1) x 1
    CW{i}(end) = 1.0;
end

models_mixed = models;

for j = 1:M
    mu_ij = Pi(:, j) .* mu ./ c_bar(j);
    mu_ij = max(mu_ij, 0);
    mu_ij = mu_ij / sum(mu_ij);
    cmu   = [0; cumsum(mu_ij(:))];
    cmu(end) = 1.0;

    new_p = zeros(nx, N);
    for m = 1:N
        % 选源模型
        r1    = rand();
        i_src = find(cmu(2:end) >= r1, 1, 'first');
        if isempty(i_src); i_src = M; end

        % 从源模型粒子集重采样
        r2    = rand();
        l_src = find(CW{i_src}(2:end) >= r2, 1, 'first');
        if isempty(l_src); l_src = N; end

        new_p(:, m) = models(i_src).particles(:, l_src);
    end

    models_mixed(j).particles = new_p;
    models_mixed(j).weights   = ones(1, N) / N;
    models_mixed(j) = pf_extract_state(models_mixed(j));
end
end
