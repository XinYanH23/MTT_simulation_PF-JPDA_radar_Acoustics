function [track, frame_stats] = track_imm_pda_update_ab( ...
    track, Z, R_list, assoc, t_idx, cfg, ac_pack, use_prod_pf, ac_field_mode)
%TRACK_IMM_PDA_UPDATE_AB  A/B 验证专用轨迹更新。
%   use_prod_pf=true  → 调用生产 pf_jpda_update（B: L_softmax）
%   use_prod_pf=false → 调用 pf_jpda_update_ab（A: post 等）

if nargin < 8 || isempty(use_prod_pf)
    use_prod_pf = false;
end
if nargin < 9
    ac_field_mode = 'post';
end

H    = cfg.H;
M    = length(track.models);

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

Lambda    = zeros(M, 1);
ess_frame = inf;
n_rs      = 0;

for j = 1:M
    if use_prod_pf
        [track.models(j), Lj, ess_j] = pf_jpda_update( ...
            track.models(j), H, ...
            Z_valid, R_valid, beta_valid, beta0_val, cfg.lambda_c, ...
            ac_pack, cfg);
        n_rs = n_rs + (ess_j < cfg.pf_resample_thresh * track.models(j).N);
    else
        [track.models(j), Lj, ess_j, dbg] = pf_jpda_update_ab( ...
            track.models(j), H, ...
            Z_valid, R_valid, beta_valid, beta0_val, cfg.lambda_c, ...
            ac_pack, cfg, ac_field_mode);
        n_rs = n_rs + dbg.resampled;
    end
    Lambda(j) = Lj;
    ess_frame = min(ess_frame, ess_j);
end

if isfield(track, 'c_bar_')
    c_bar = track.c_bar_;
else
    c_bar = track.Pi' * track.mu;
end
c_bar    = max(c_bar(:), 1e-300);
track.mu = imm_update_mu(track.mu, Lambda, c_bar);
[track.x, track.P] = imm_fuse(track.models, track.mu);

if use_pf(cfg) && isfinite(ess_frame)
    if ~isfield(track, 'ess_hist') || isempty(track.ess_hist)
        track.ess_hist = ess_frame;
    else
        track.ess_hist(end+1) = ess_frame;
    end
end

frame_stats = struct();
frame_stats.ess         = ess_frame;
frame_stats.mu          = track.mu(:)';
frame_stats.n_resample  = n_rs;
frame_stats.n_gate_meas = Dv;
frame_stats.beta0       = beta0_val;
if Dv > 0
    frame_stats.beta_d_max = max(beta_valid);
else
    frame_stats.beta_d_max = 0;
end
[~, mj] = max(track.mu);
frame_stats.dominant_model = mj;
if use_pf(cfg)
    frame_stats.particles_xy = track.models(mj).particles(1:2, :);
    frame_stats.weights      = track.models(mj).weights;
end
end

function tf = use_pf(cfg)
tf = isfield(cfg, 'use_pf') && cfg.use_pf;
end
