function assoc = jpda_run(tracks, Z, R_list, cfg)
% JPDA_RUN  完整 JPDA 流水线：门控 → 聚类 → β 计算。
%
%   assoc = jpda_run(tracks, Z, R_list, cfg)
%
%   INPUTS
%     tracks  : 1×T struct 数组（活跃轨迹）
%     Z       : 2×D  量测位置矩阵
%     R_list  : 2×2×D  每探测量测噪声
%     cfg     : config_tracker 输出
%
%   OUTPUT
%     assoc.valid_mat  : T×D logical
%     assoc.innov_data : T×D cell（新息等数据）
%     assoc.beta       : T×D  边缘关联概率
%     assoc.beta0      : T×1  漏检概率
%     assoc.clusters   : 1×K cell
%
%   调用顺序：
%     jpda_gating → jpda_cluster → jpda_compute_beta

T = length(tracks);
D = size(Z, 2);

% 默认值（无探测或无轨迹时）
assoc.valid_mat  = false(T, D);
assoc.innov_data = cell(T, D);
assoc.beta       = zeros(T, D);
assoc.beta0      = ones(T, 1);
assoc.clusters   = {};

if T == 0 || D == 0
    return
end

[assoc.valid_mat, assoc.innov_data] = jpda_gating(tracks, Z, R_list, cfg);
assoc.clusters = jpda_cluster(assoc.valid_mat);
[assoc.beta, assoc.beta0] = jpda_compute_beta( ...
    assoc.valid_mat, assoc.innov_data, assoc.clusters, cfg);
end
