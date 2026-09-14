function pack = build_acoustic_prob_field(measured_lp, cfg, method, prev_post)
%BUILD_ACOUSTIC_PROB_FIELD 声学概率场量测主入口（供后续 MATLAB 跟踪器引入）。
%
%   pack = BUILD_ACOUSTIC_PROB_FIELD(measured_lp, cfg)
%   pack = BUILD_ACOUSTIC_PROB_FIELD(measured_lp, cfg, method, prev_post)
%
%   给定各节点实测/仿真 SPL，输出该帧声学概率场量测包。这是声学侧唯一对外接口；
%   不输出坐标量测 z=x+noise，仅输出概率场与点估计供数据关联/融合层查询。
%
%   输入：
%     measured_lp : Mx1 各节点 SPL (dB)，缺测填 NaN
%     cfg         : acoustic_config()
%     method      : 'softmax' | 'intensity' | 'nms_gmm'（默认 'softmax'）
%                   决定送入递归贝叶斯的似然场（三种场始终全部计算并返回）
%     prev_post   : 上一帧后验（递归贝叶斯用）；首帧传 []
%
%   输出 struct pack 字段：
%     xc, zc, XX, ZZ          网格坐标 (mm)
%     S                        融合对数似然得分场
%     L_rel                    相对似然场 exp(S-max S)，论文主线查询对象
%     L_softmax,L_intensity,L_nms_gmm   三种工程映射（非主理论）
%     peaks                    NMS 峰值 [x,z,score,weight]
%     post                     递归贝叶斯后验（基于 method 选定的似然场）
%     bayes_info               熵/是否更新
%     est_mean, est_map        后验点估计 [x_mm,z_mm]
%     snr_db, weights          各节点 SNR 与可靠性权重
%     p_detect_m, p_detect_joint, detected   探测概率与判定

if nargin < 3 || isempty(method); method = 'softmax'; end
if nargin < 4; prev_post = []; end

% 1. 节点 SNR / 可靠性权重 / 探测概率
[snr_db, w, p_detect_m, p_detect_joint] = node_snr_weights(measured_lp, cfg);
valid = ~isnan(snr_db(:)) & (w(:) > 0);
if any(valid)
    snr_bar = mean(snr_db(valid));
else
    snr_bar = -Inf;
end

% 2. 多节点协同对数似然得分场
[S, xc, zc, XX, ZZ] = build_score_grid(measured_lp, snr_db, w, cfg);

% 3. 相对似然场（论文主线）与三种工程映射（仅保留供可视化/旧实验）
Smax = max(S(:));
if isfinite(Smax)
    L_rel = exp(S - Smax);
else
    L_rel = ones(size(S));
end
L_softmax   = prob_field_softmax(S, snr_bar, cfg);
L_intensity = prob_field_intensity(S, cfg);
[L_nms_gmm, peaks] = prob_field_nms_gmm(S, XX, ZZ, snr_bar, cfg);

% 4. 选定似然场做递归贝叶斯
switch lower(method)
    case 'softmax';   L_sel = L_softmax;
    case 'intensity'; L_sel = L_intensity;
    case 'nms_gmm';   L_sel = L_nms_gmm;
    otherwise; error('build_acoustic_prob_field:method', '未知 method: %s', method);
end
[post, bayes_info] = recursive_bayes_field(prev_post, L_sel, cfg);

% 5. 点估计
[est_mean, est_map] = field_point_estimate(post, XX, ZZ);

pack = struct( ...
    'xc', xc, 'zc', zc, 'XX', XX, 'ZZ', ZZ, ...
    'S', S, 'L_rel', L_rel, ...
    'L_softmax', L_softmax, 'L_intensity', L_intensity, 'L_nms_gmm', L_nms_gmm, ...
    'peaks', peaks, ...
    'post', post, 'bayes_info', bayes_info, ...
    'est_mean', est_mean, 'est_map', est_map, ...
    'snr_db', snr_db, 'weights', w, ...
    'p_detect_m', p_detect_m, 'p_detect_joint', p_detect_joint, ...
    'detected', p_detect_joint >= cfg.det_threshold, ...
    'method', method, 'snr_bar_db', snr_bar);

end
