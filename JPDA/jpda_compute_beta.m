function [beta, beta0] = jpda_compute_beta(valid_mat, innov_data, clusters, cfg)
% JPDA_COMPUTE_BETA  计算 JPDA 边缘关联概率 β_{t,d} 和漏检概率 β_{t,0}。
%
%   [beta, beta0] = jpda_compute_beta(valid_mat, innov_data, clusters, cfg)
%
%   INPUTS
%     valid_mat  : T×D logical        验证矩阵
%     innov_data : T×D cell           新息数据（来自 jpda_gating）
%     clusters   : 1×K cell           簇列表（来自 jpda_cluster）
%     cfg        : config_tracker     含 Pd, lambda_c, max_events
%
%   OUTPUTS
%     beta  : T×D double  β_{t,d}  目标 t 来自探测 d 的边缘概率
%     beta0 : T×1 double  β_{t,0}  目标 t 本帧漏检的边缘概率
%
%   ── 算法（Bar-Shalom 2001, §6.2.2 ~ 6.2.3）──────────────────────
%
%   联合事件概率（未归一化对数形式）：
%
%     log P(θ) ∝  Σ_{t:θ(t)≠0} [ log Pd + log L_{t,θ(t)} − log λ_c ]
%              +  Σ_{t:θ(t)=0} log(1 − Pd)
%
%   其中 L_{t,d} = N(v_{t,d}; 0, S_{t,d}) 为高斯量测似然。
%
%   归一化后，边缘化得到：
%
%     β_{t,d} = Σ_{θ: θ(t)=d} P(θ)   / Σ_θ P(θ)
%     β_{t,0} = Σ_{θ: θ(t)=0} P(θ)   / Σ_θ P(θ)
%
%   满足 β_{t,0} + Σ_d β_{t,d} = 1，∀ t。
%
%   注意：与"单目标 PDA"不同，真正的 JPDA 在联合事件层面强制
%   "每个探测最多来自一个目标"的约束，使得不同目标的 β 在
%   分母中相互耦合，从而正确处理密集目标场景。

[T, D] = size(valid_mat);

Pd          = cfg.Pd;
lc          = cfg.lambda_c;
log_Pd_lc   = log(Pd) - log(lc);
log_1_Pd    = log(1 - Pd + 1e-300);

beta  = zeros(T, D);
beta0 = zeros(T, 1);

%% ── 预计算高斯量测对数似然 ────────────────────────────────────────
% logL(t,d) = log N(v;0,S) = −0.5 v'S⁻¹v − 0.5 log|2πS|
logL = -inf(T, D);
for t = 1 : T
    for d = 1 : D
        if valid_mat(t, d)
            nd        = innov_data{t, d};
            logL(t,d) = -0.5 * nd.d2 - 0.5 * log(det(2*pi * nd.S) + 1e-300);
        end
    end
end

%% ── 逐簇枚举与边缘化 ─────────────────────────────────────────────
for k = 1 : length(clusters)
    t_idx = clusters{k}.targets;
    d_idx = clusters{k}.detections;
    Tloc  = length(t_idx);
    Dloc  = length(d_idx);

    % 无探测的簇：所有目标漏检
    if Dloc == 0
        for ti = 1 : Tloc
            beta0(t_idx(ti)) = 1.0;
        end
        continue
    end

    % 局部验证子矩阵
    valid_loc = valid_mat(t_idx, d_idx);
    logL_loc  = logL(t_idx, d_idx);       % Tloc×Dloc

    % 枚举可行联合事件
    Theta    = jpda_enumerate(Tloc, Dloc, valid_loc, cfg.max_events);
    n_events = length(Theta);

    % 计算每个事件的未归一化对数概率
    log_P = zeros(n_events, 1);
    for e = 1 : n_events
        theta = Theta{e};    % 1×Tloc，值域 0..Dloc
        lp    = 0;
        for ti = 1 : Tloc
            if theta(ti) == 0
                lp = lp + log_1_Pd;
            else
                lp = lp + log_Pd_lc + logL_loc(ti, theta(ti));
            end
        end
        log_P(e) = lp;
    end

    % 数值稳定归一化（减去最大值后 exp）
    log_P = log_P - max(log_P);
    P_evt = exp(log_P);
    P_evt = P_evt / sum(P_evt);

    % 边缘化
    for e = 1 : n_events
        theta = Theta{e};
        pw    = P_evt(e);
        for ti = 1 : Tloc
            t = t_idx(ti);
            if theta(ti) == 0
                beta0(t) = beta0(t) + pw;
            else
                d          = d_idx(theta(ti));
                beta(t, d) = beta(t, d) + pw;
            end
        end
    end
end

%% ── 补全孤立目标（无任何门内探测）的漏检概率 ──────────────────────
for t = 1 : T
    if ~any(valid_mat(t, :)) && beta0(t) < 1e-12
        beta0(t) = 1.0;
    end
end
end
