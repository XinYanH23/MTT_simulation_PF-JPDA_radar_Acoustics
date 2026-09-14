function R_bar = compute_R_bar(beta_valid, R_valid, beta0, normalize_Rbar)
%COMPUTE_R_BAR  PDA 合成量测噪声（答辩 B2 对照入口）。
%
%   normalize_Rbar = false : R̄ = Σ_d β_d R_d
%   normalize_Rbar = true  : R̄ = Σ_d β_d R_d / (1 − β_0)
%
%   当 β_0 ≈ 0 时两者几乎相同；高漏检/多杂波时用于实验区分。

    beta_valid = beta_valid(:)';
    Dv         = size(R_valid, 3);
    m          = size(R_valid, 1);
    R_bar      = zeros(m, m);

    if Dv < 1
        return
    end

    beta_sum = sum(beta_valid);
    if beta_sum < 1e-12
        R_bar = R_valid(:, :, 1);
        return
    end

    for d = 1 : Dv
        R_bar = R_bar + beta_valid(d) * R_valid(:, :, d);
    end

    if normalize_Rbar
        denom = max(1 - beta0, 1e-6);
        R_bar = R_bar / denom;
    end
end
