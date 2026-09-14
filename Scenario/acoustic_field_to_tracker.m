function [ac_pack, prev_post_out] = acoustic_field_to_tracker(src_pos_m, ac_cfg, prev_post, method)
%ACOUSTIC_FIELD_TO_TRACKER  声学概率场生成接口（主循环每帧调用）。
%
%   多目标策略（与论文“空间似然场可多峰”一致）：
%     对每个声源单独建场 L_i，再取 L = max_i L_i（峰值保留、避免 Softmax
%     全局归一化把次强峰压到数值零）。联合探测概率取各源 pdj 的并。
%
%   坐标系：
%     输入 src_pos_m : [x_east, y_north, z_height] (m)
%     声学 mic/src   : [X_east, Y_height, Z_north] (mm)

if nargin < 3;  prev_post = [];    end
if nargin < 4;  method    = 'softmax';  end

if isempty(src_pos_m)
    ac_pack      = struct('detected', false, 'post', [], ...
                          'p_detect_joint', 0, 'method', method);
    prev_post_out = prev_post;
    return
end

% 多目标场景下提高 Softmax 温度，避免峰过尖
if ~isfield(ac_cfg, 'softmax_tau0_set')
    ac_cfg.softmax_tau0 = max(ac_cfg.softmax_tau0, 3.0);
    ac_cfg.softmax_tau_c = 0.05;
end

K = size(src_pos_m, 1);
packs = cell(K, 1);
pdj_all = zeros(K, 1);

for i = 1:K
    src_xyz_mm = src_to_ac_mm(src_pos_m(i, :), ac_cfg);
    M = size(ac_cfg.mic_coords, 1);
    measured_lp = zeros(M, 1);
    for m = 1:M
        measured_lp(m) = single_mic_spl(src_xyz_mm, m, ac_cfg);
    end
    if isfield(ac_cfg, 'sigma0_db') && ac_cfg.sigma0_db > 0
        measured_lp = measured_lp + 0.3 * ac_cfg.sigma0_db * randn(M, 1);
    end
    % 单源场：不做跨帧贝叶斯（多源混合后再传 prev_post 易串扰）
    packs{i} = build_acoustic_prob_field(measured_lp, ac_cfg, method, []);
    pdj_all(i) = packs{i}.p_detect_joint;
end

% 以最强源的 pack 为模板，工程场取逐点 max（理论主线只用单帧 L_rel）
ac_pack = packs{1};
L_names = {'L_rel', 'L_softmax', 'L_intensity', 'L_nms_gmm'};
for ni = 1:numel(L_names)
    nm = L_names{ni};
    if ~isfield(ac_pack, nm) || isempty(ac_pack.(nm)), continue; end
    Lmax = ac_pack.(nm);
    for i = 2:K
        if isfield(packs{i}, nm) && ~isempty(packs{i}.(nm)) ...
                && isequal(size(packs{i}.(nm)), size(Lmax))
            Lmax = max(Lmax, packs{i}.(nm));
        end
    end
    ac_pack.(nm) = Lmax;
end

% 联合探测：至少一源被探测
ac_pack.p_detect_joint = 1 - prod(1 - pdj_all);
ac_pack.detected = ac_pack.p_detect_joint >= ac_cfg.det_threshold;
ac_pack.method = method;
ac_pack.n_sources = K;

% 递归后验：用混合后的主似然做一步贝叶斯（可选）
if ~isempty(prev_post) && isfield(ac_pack, 'L_softmax') && ~isempty(ac_pack.L_softmax)
    try
        ac_pack.post = recursive_bayes_field(prev_post, ac_pack.L_softmax, ac_cfg);
    catch
        ac_pack.post = ac_pack.L_softmax;
    end
else
    ac_pack.post = ac_pack.L_softmax;
end
prev_post_out = ac_pack.post;
end

% -------------------------------------------------------------------------
function xyz = src_to_ac_mm(src_pos_m, ac_cfg)
X_mm = src_pos_m(1) * 1000;
Z_mm = src_pos_m(2) * 1000;
if numel(src_pos_m) >= 3 && isfinite(src_pos_m(3)) && src_pos_m(3) > 1
    Y_mm = src_pos_m(3) * 1000;
else
    Y_mm = ac_cfg.src_height_mm;
end
xyz = [X_mm, Y_mm, Z_mm];
end
