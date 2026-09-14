function track = track_imm_pf_update(track, Z, R_list, assoc, t_idx, ac_pack, cfg)
%TRACK_IMM_PF_UPDATE  单条轨迹 IMM-PF + JPDA 后验混合 + 声学似然修正。
%
%   track = TRACK_IMM_PF_UPDATE(track, Z, R_list, assoc, t_idx, ac_pack, cfg)
%
%   调用约定：调用前主循环已完成 imm_pf_mix + imm_pf_predict，
%             track.models 已是预测后粒子集，track.c_bar_ 已暂存。
%
%   INPUTS
%     track   : 轨迹结构体（含 .models[j].particles/.weights）
%     Z       : 2×D    全部量测
%     R_list  : 2×2×D  全部量测噪声
%     assoc   : jpda_run 输出（beta, beta0, valid_mat）—— 不修改
%     t_idx   : 本轨迹在 assoc 矩阵中的行索引
%     ac_pack : 声学概率场 struct（build_acoustic_prob_field 输出），空 [] 则跳过
%     cfg     : config_tracker 输出（含 pf_* 字段）
%
%   OUTPUTS
%     track : 更新后（模型粒子权重已更新，mu 已更新，x/P 已融合）

H  = cfg.H;
M  = length(track.models);
D  = size(Z, 2);

% ── 提取本轨迹 JPDA 关联结果 ──────────────────────────────────────────────
beta_row  = assoc.beta(t_idx, :);
beta0_val = assoc.beta0(t_idx);
d_valid   = find(assoc.valid_mat(t_idx, :));
Dv        = length(d_valid);

Z_valid    = zeros(2, Dv);
R_valid    = zeros(2, 2, Dv);
beta_valid = zeros(1, Dv);
for k = 1:Dv
    d              = d_valid(k);
    Z_valid(:, k)  = Z(:, d);
    R_valid(:,:,k) = R_list(:, :, d);
    beta_valid(k)  = beta_row(d);
end

% ── 对每个模型执行 PF 权重更新 ────────────────────────────────────────────
Lambda = zeros(M, 1);

for j = 1:M
    [track.models(j), Lambda(j)] = pf_jpda_update( ...
        track.models(j), H, ...
        Z_valid, R_valid, beta_valid, beta0_val, ...
        cfg.lambda_c, ac_pack, cfg);
end

% ── IMM 模型概率更新（接口不变）──────────────────────────────────────────
if isfield(track, 'c_bar_')
    c_bar = track.c_bar_;
else
    c_bar = track.Pi' * track.mu;
end
c_bar    = max(c_bar(:), 1e-300);
track.mu = imm_update_mu(track.mu, Lambda, c_bar);

% ── IMM 融合（接口不变，使用粒子均值/协方差）──────────────────────────────
[track.x, track.P] = imm_fuse(track.models, track.mu);
end
