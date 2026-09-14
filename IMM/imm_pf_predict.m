function models = imm_pf_predict(models)
%IMM_PF_PREDICT  粒子传播（替代 KF 版 imm_predict）。
%
%   models = IMM_PF_PREDICT(models)
%
%   每个模型独立传播：x^(m) = F*x^(m) + noise,  noise ~ N(0,Q)
%   Cholesky 一次采样全部粒子（向量化，无逐粒子循环）。
%   传播后同步 .x .P（供 imm_fuse 门控使用）。

for j = 1:length(models)
    F = models(j).F;
    Q = models(j).Q;
    N = models(j).N;

    L     = chol(Q, 'lower');
    noise = L * randn(size(Q, 1), N);

    models(j).particles = F * models(j).particles + noise;
    models(j) = pf_extract_state(models(j));
end
end
