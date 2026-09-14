function [x_upd, P_upd, Lambda] = ekf_pda_update(x_pred, P_pred, H, ...
                                                   Z_valid, R_valid, ...
                                                   beta_valid, beta0, lambda_c, ...
                                                   normalize_Rbar)
% EKF_PDA_UPDATE  完整 PDA 状态/协方差更新（含 innovation-spread 项）。
%
%   [x_upd, P_upd, Lambda] = ekf_pda_update(x_pred, P_pred, H, ...
%                                Z_valid, R_valid, beta_valid, beta0, lambda_c)
%
%   INPUTS
%     x_pred      : n×1   预测状态
%     P_pred      : n×n   预测协方差
%     H           : m×n   量测矩阵
%     Z_valid     : m×Dv  门内探测（Dv 列）
%     R_valid     : m×m×Dv  各探测量测噪声
%     beta_valid  : 1×Dv  β_{t,d}（门内探测的关联概率）
%     beta0       : scalar  β_{t,0}（漏检概率）
%     lambda_c    : scalar  杂波空间密度（用于模型似然计算）
%
%   OUTPUTS
%     x_upd  : n×1  更新后状态
%     P_upd  : n×n  更新后协方差（完整 PDA 三项式）
%     Lambda : scalar  模型似然（供 imm_update_mu 使用）
%
%   ── 完整 PDA 公式（Bar-Shalom 2001, §6.3.2）──────────────────────
%
%   定义：
%     v_d  = z_d − H x̂(−)           （每个门内探测的新息）
%     v̄   = Σ_d β_d v_d             （加权合成新息）
%     R̄   = Σ_d β_d R_d             （默认，Bar-Shalom §6.3.2）
%     或   R̄ = Σ_d β_d R_d / (1−β_0) （normalize_Rbar=true，答辩 B2 对照）
%     S̄   = H P(−) H' + R̄
%     K    = P(−) H' S̄⁻¹           （卡尔曼增益）
%
%   状态更新：
%     x̂(+) = x̂(−) + K v̄
%
%   协方差更新（三项）：
%     项 1（漏检分量）：         β_0 · P(−)
%     项 2（条件协方差分量）：   (1 − β_0) · P_c,    P_c = (I − KH) P(−)
%     项 3（innovation spread）： K · [Σ_d β_d v_d v_d' − v̄ v̄'] · K'
%
%     P(+) = β_0 P(−) + (1−β_0) P_c + K [Σ_d β_d v_d v_d' − v̄ v̄'] K'
%
%   ────────────────────────────────────────────────────────────────────
%   ★ 物理意义（必须理解，否则容易误删）：
%
%   项 1：若目标被漏检（β_0 > 0），此部分保留预测协方差，
%         体现"没有量测，不确定性不减小"的物理事实。
%
%   项 2：若目标被某探测正确检测，则条件协方差 P_c < P(−)，
%         体现量测带来的信息增益。(1−β_0) 是检测概率的平均。
%
%   项 3（最容易被遗漏！）：
%         当存在多个候选探测、真实来源不确定时，
%         状态估计的方差需额外增加"因不知道用哪个探测更新"
%         所引入的不确定性。若 β 集中于一个探测，此项趋近于零；
%         若 β 均匀分散在多个探测，此项最大。
%         省略此项会使滤波器在密集环境下协方差虚假偏小，
%         导致后续帧门控过窄、轨迹丢失。
%   ────────────────────────────────────────────────────────────────────

if nargin < 9 || isempty(normalize_Rbar)
    normalize_Rbar = false;
end

Dv = size(Z_valid, 2);
n  = size(x_pred, 1);
m  = size(H, 1);

%% ── 无门内探测：仅做预测传播（漏检帧）──────────────────────────────
if Dv == 0
    x_upd  = x_pred;
    P_upd  = P_pred;
    Lambda = lambda_c;   % 漏检时似然等于杂波密度（最低基准值）
    return
end

%% ── β 加权量测噪声 R̄（B2：两种实现可切换）────────────────────────
beta_valid = beta_valid(:)';        % 强制行向量 1×Dv
R_bar      = compute_R_bar(beta_valid, R_valid, beta0, normalize_Rbar);

%% ── 新息协方差与卡尔曼增益 ──────────────────────────────────────────
S_bar = H * P_pred * H' + R_bar;
S_bar = (S_bar + S_bar') / 2;
K     = P_pred * H' / S_bar;

%% ── 新息矩阵 V（每列对应一个门内探测）───────────────────────────────
V = zeros(m, Dv);
for d = 1 : Dv
    V(:, d) = Z_valid(:, d) - H * x_pred;
end

%% ── 合成新息 v̄ ───────────────────────────────────────────────────
v_bar = V * beta_valid(:);      % m×1

%% ── 状态更新 ────────────────────────────────────────────────────────
x_upd = x_pred + K * v_bar;

%% ── 协方差更新（三项 PDA 公式）──────────────────────────────────────

% 项 1：漏检分量
P_term1 = beta0 * P_pred;

% 项 2：条件协方差分量
P_c     = P_pred - K * S_bar * K';
P_c     = (P_c + P_c') / 2;
P_term2 = (1 - beta0) * P_c;

% 项 3：innovation spread（不可省略）
S_spread = zeros(m, m);
for d = 1 : Dv
    S_spread = S_spread + beta_valid(d) * (V(:, d) * V(:, d)');
end
P_term3 = K * (S_spread - v_bar * v_bar') * K';

P_upd = P_term1 + P_term2 + P_term3;
P_upd = (P_upd + P_upd') / 2;

%% ── 模型似然 Λ（供 IMM 模型概率更新）──────────────────────────────
% Λ ≈ λ_c + Pd · Σ_d L_d
% 其中 L_d = N(v_d; 0, S̄)（用 S_bar 近似，避免重复计算）
%
% 此值仅用于 imm_update_mu 中各模型之间的相对比较，
% 绝对值不影响结果，但应随门内探测数量单调增大。
L_sum = 0;
inv_S_bar = S_bar \ eye(m);    % 一次计算，多次复用
det_S_bar = det(2 * pi * S_bar);
for d = 1 : Dv
    v   = V(:, d);
    L_d = exp(-0.5 * v' * inv_S_bar * v) / sqrt(max(det_S_bar, 1e-300));
    L_sum = L_sum + beta_valid(d) * L_d;
end

Lambda = max(lambda_c + L_sum, 1e-300);
end
