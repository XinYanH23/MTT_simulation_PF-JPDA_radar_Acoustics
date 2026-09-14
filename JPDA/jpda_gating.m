function [valid_mat, innov_data] = jpda_gating(tracks, Z, R_list, cfg)
% JPDA_GATING  马氏距离验证门（Validation Gate）。
%
%   [valid_mat, innov_data] = jpda_gating(tracks, Z, R_list, cfg)
%
%   INPUTS
%     tracks   : 1×T struct 数组（活跃轨迹，含 .x .P，均为预测后融合值）
%     Z        : 2×D  量测矩阵（每列一个探测的位置量测）
%     R_list   : 2×2×D  每个探测的量测噪声协方差
%     cfg      : config_tracker 输出
%
%   OUTPUTS
%     valid_mat  : T×D logical  true 表示 (t,d) 通过验证门
%     innov_data : T×D cell     每个通过门的 (t,d) 存储：
%                    .v      新息向量 (2×1)
%                    .S      新息协方差 (2×2)
%                    .S_inv  S 的逆 (2×2)
%                    .K      卡尔曼增益 (4×2)，用融合 P 计算
%                    .d2     马氏距离²
%
%   门限条件（Bar-Shalom §3.4.1）：
%     d² = v' S⁻¹ v < γ
%   其中 v = z_d − H x̂_t(−)，S = H P_t(−) H' + R_d，γ = χ²(n_z, p)

H   = cfg.H;
g   = cfg.gate_chi2;
T   = length(tracks);
D   = size(Z, 2);

valid_mat  = false(T, D);
innov_data = cell(T, D);

for t = 1 : T
    x_pred = tracks(t).x;        % 4×1 融合预测状态
    P_pred = tracks(t).P;        % 4×4 融合预测协方差

    for d = 1 : D
        R = R_list(:, :, d);
        z = Z(:, d);

        v = z - H * x_pred;
        S = H * P_pred * H' + R;
        S = (S + S') / 2;

        % 使用 Cholesky 分解计算马氏距离，数值更稳定
        [L_chol, flag] = chol(S, 'lower');
        if flag ~= 0
            % 矩阵非正定，加正则化
            S    = S + 1e-6 * eye(size(S));
            S_inv = inv(S);
        else
            S_inv = (L_chol' \ (L_chol \ eye(size(S))));
        end

        d2 = v' * S_inv * v;

        if isfinite(d2) && d2 < g
            valid_mat(t, d)  = true;
            K                = P_pred * H' * S_inv;
            innov_data{t, d} = struct('v', v, 'S', S, 'S_inv', S_inv, ...
                                      'K', K, 'd2', d2);
        end
    end
end
end
