function track = track_imm_pda_update(track, Z, R_list, assoc, t_idx, cfg, ac_pack)
% NOTE: 调用前主循环已完成 imm_mix + imm_predict，
%       track.models 已是预测后状态，track.c_bar_ 已由主循环暂存。
% TRACK_IMM_PDA_UPDATE  对单条轨迹执行完整的 IMM + PDA 更新。
%
%   track = track_imm_pda_update(track, Z, R_list, assoc, t_idx, cfg)
%
%   INPUTS
%     track   : 轨迹结构体（已完成 IMM 混合与预测）
%     Z       : 2×D  全部量测
%     R_list  : 2×2×D  全部量测噪声
%     assoc   : jpda_run 输出（含 beta, beta0, valid_mat, innov_data）
%     t_idx   : 本轨迹在 assoc 矩阵中的行索引
%     cfg     : config_tracker 输出
%
%   OUTPUT
%     track : 更新后的轨迹（models, mu, x, P 均已更新）
%
%   流程：
%     1. 从 assoc 中提取本轨迹的 β 和门内探测
%     2. 对每个 IMM 模型：PF 走 pf_jpda_update；KF 走 ekf_pda_update
%        后接 ekf_acoustic_field_update（同一 L_k、同一 γ 混合）
%     3. 调用 imm_update_mu 更新模型概率
%     4. 调用 imm_fuse 得到融合估计

if nargin < 7;  ac_pack = [];  end
use_pf = isfield(cfg, 'use_pf') && cfg.use_pf;

H    = cfg.H;
M    = length(track.models);
D    = size(Z, 2);

% 本轨迹的 β 行
beta_row  = assoc.beta(t_idx, :);          % 1×D
beta0_val = assoc.beta0(t_idx);            % scalar

% 找出门内探测（valid_mat 中本行为 true 的列）
d_valid = find(assoc.valid_mat(t_idx, :)); % 1×Dv

% 提取门内探测的量测与噪声
Dv      = length(d_valid);
Z_valid = zeros(2, Dv);
R_valid = zeros(2, 2, Dv);
beta_valid = zeros(1, Dv);

for k = 1 : Dv
    d              = d_valid(k);
    Z_valid(:, k)  = Z(:, d);
    R_valid(:,:,k) = R_list(:, :, d);
    beta_valid(k)  = beta_row(d);
end

%% ── 对每个模型执行 PDA 更新（KF 或 PF 分支） ────────────────────
Lambda = zeros(M, 1);
ess_frame = inf;

for j = 1 : M
    if use_pf
        %% ── PF 分支 ────────────────────────────────────────────
        [track.models(j), Lj, ess_j] = pf_jpda_update( ...
            track.models(j), H, ...
            Z_valid, R_valid, beta_valid, beta0_val, cfg.lambda_c, ...
            ac_pack, cfg);
        Lambda(j) = Lj;
        ess_frame = min(ess_frame, ess_j);

    else
        %% ── KF 分支：雷达 PDA + 与 PF 相同的声学场乘性更新 ─────
        x_j = track.models(j).x;
        P_j = track.models(j).P;

        if Dv == 0
            x_upd = x_j;
            P_upd = P_j;
            Lj    = cfg.lambda_c;
        else
            R_valid_j = R_valid;
            norm_R    = isfield(cfg, 'pda_Rbar_normalized') && cfg.pda_Rbar_normalized;
            [x_upd, P_upd, Lj] = ekf_pda_update( ...
                x_j, P_j, H, Z_valid, R_valid_j, beta_valid, beta0_val, ...
                cfg.lambda_c, norm_R);
        end

        % 与 PF 相同：ℓ_fused = ℓ_radar · [γ L(x)+(1-γ)]
        % 高斯无法逐粒子积分多峰场，用 σ 点再加权匹配一、二阶矩。
        if ~isempty(ac_pack)
            [x_upd, P_upd, ell_ac_bar] = ekf_acoustic_field_update( ...
                x_upd, P_upd, ac_pack, cfg);
            Lj = Lj * ell_ac_bar;
        end
        track.models(j).x = x_upd;
        track.models(j).P = P_upd;
        Lambda(j)         = Lj;
    end
end

%% ── IMM 模型概率更新 ─────────────────────────────────────────────
% c_bar 由主循环在 imm_mix 后暂存于 track.c_bar_
% （c̄_j = Σ_i π_{ij} μ_i，在混合步骤确定，此后不变）
if isfield(track, 'c_bar_')
    c_bar = track.c_bar_;
else
    c_bar = track.Pi' * track.mu;   % 回退：不应执行到此
end
c_bar    = max(c_bar(:), 1e-300);
track.mu = imm_update_mu(track.mu, Lambda, c_bar);

%% ── IMM 融合输出 ─────────────────────────────────────────────────
[track.x, track.P] = imm_fuse(track.models, track.mu);

if use_pf && isfinite(ess_frame)
    if ~isfield(track, 'ess_hist') || isempty(track.ess_hist)
        track.ess_hist = ess_frame;
    else
        track.ess_hist(end+1) = ess_frame; %#ok<AGROW>
    end
end
end
