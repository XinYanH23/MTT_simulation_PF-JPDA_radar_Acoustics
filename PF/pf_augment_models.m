function models = pf_augment_models(models, cfg)
%PF_AUGMENT_MODELS  为已有 IMM models 添加粒子字段（不修改 x/P/F/Q）。
%
%   models = PF_AUGMENT_MODELS(models, cfg)
%
%   调用时机：track_create 之后，首次 imm_pf_mix 之前。
%   每个模型的粒子从当前 model.x / model.P 的高斯分布初始化。

N = cfg.pf_N;

for j = 1:length(models)
    pf_j = pf_init(models(j).x, models(j).P, N, cfg.pf_init_cov_scale);
    models(j).particles = pf_j.particles;
    models(j).weights   = pf_j.weights;
    models(j).N         = pf_j.N;
end
end
