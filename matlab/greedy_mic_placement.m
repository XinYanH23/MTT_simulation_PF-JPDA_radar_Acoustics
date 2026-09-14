function result = greedy_mic_placement(candidates, cfg, cand_A, cand_n, cand_nf)
%GREEDY_MIC_PLACEMENT 贪婪节点布设算法：以全域平均探测置信度为目标函数，
%   在覆盖率满足要求的前提下选出节点数最少的最优布局。
%
%   result = GREEDY_MIC_PLACEMENT(candidates, cfg)
%   result = GREEDY_MIC_PLACEMENT(candidates, cfg, cand_A, cand_n, cand_nf)
%
%   输入：
%     candidates : Kx3 候选节点世界坐标 [X,Y,Z] (mm)
%     cfg        : acoustic_config() 返回的参数 struct（含 §8 布局优化参数）
%     cand_A     : Kx1 各候选节点 A 参数；省略则用 cfg.placement_A_default
%     cand_n     : Kx1 各候选节点 n 参数；省略则用 cfg.placement_n_default
%     cand_nf    : Kx1 各候选节点环境底噪 (dB)；省略则用 cfg.placement_nf_default
%
%   算法（贪婪）：
%     初始：已选集 S = {}，剩余候选集 R = {1..K}
%     每步：
%       对所有 c ∈ R，计算 S∪{c} 的全域平均 P_D_joint
%       选 c* = argmax mean_pd(S∪{c})
%       S ← S∪{c*}，R ← R\{c*}
%     停止：mean_pd(S) >= min_coverage  OR  |S| >= max_nodes
%
%   目标函数：mean_pd = mean_{u} P_D_joint(u)
%   覆盖率指标：coverage_frac = E[P_D_joint(u) >= threshold]
%
%   输出 struct result：
%     sel_idx        : 1xM 已选节点在 candidates 中的行号（选取顺序）
%     sel_coords     : Mx3 已选节点坐标
%     mean_pd_hist   : 1x(M+1) 每步加入后的全域平均 P_D_joint（含空集=0）
%     coverage_hist  : 1x(M+1) 每步 coverage_frac（含空集=0）
%     pd_map_final   : nx x nz 最终联合探测概率场
%     marginal_gain  : 1xM 每步加入带来的 mean_pd 增量（边际贡献）
%     xc, zc         : 覆盖评估网格中心 (mm)
%     n_nodes        : 最终选取节点数
%     satisfied      : 是否在预算内达到覆盖目标

K = size(candidates, 1);

if nargin < 3 || isempty(cand_A)
    cand_A  = repmat(cfg.placement_A_default,  K, 1);
end
if nargin < 4 || isempty(cand_n)
    cand_n  = repmat(cfg.placement_n_default,  K, 1);
end
if nargin < 5 || isempty(cand_nf)
    cand_nf = repmat(cfg.placement_nf_default, K, 1);
end
cand_A  = cand_A(:);
cand_n  = cand_n(:);
cand_nf = cand_nf(:);

% --- 覆盖评估网格（计算一次，共享）---
[xc, zc]  = placement_grid_local(cfg);
[XX, ZZ]  = ndgrid(xc, zc);
grid_size = numel(XX);

% --- 预计算每个候选节点在整个网格上的 P_D 场（节省重复计算）---
PD_all = zeros(grid_size, K);          % 每列 = 候选 k 在网格上的 P_D_m
for k = 1:K
    Lp_k = lp_grid_local(candidates(k,:), cand_A(k), cand_n(k), XX, ZZ, cfg);
    snr_k = Lp_k - cand_nf(k);
    PD_all(:, k) = acoustic_detection_prob(snr_k(:), cfg);
end

% --- 贪婪主循环 ---
selected  = zeros(1, 0);       % 已选 idx
remaining = 1:K;                % 剩余候选
% 当前已选集的 (1 - P_D_m) 乘积，维护以避免每步重算
one_minus_cur = ones(grid_size, 1);

mean_pd_hist   = zeros(1, cfg.placement_max_nodes + 1);
coverage_hist  = zeros(1, cfg.placement_max_nodes + 1);
marginal_gain  = zeros(1, cfg.placement_max_nodes);

% 空集统计
mean_pd_hist(1)  = 0;
coverage_hist(1) = 0;

satisfied = false;
for step = 1:cfg.placement_max_nodes

    if isempty(remaining); break; end

    % 对每个剩余候选计算加入后的 mean_pd（向量化）
    one_minus_try = one_minus_cur .* (1 - PD_all(:, remaining));   % grid_size x |remaining|
    pd_joint_try  = 1 - one_minus_try;                              % grid_size x |remaining|
    mean_pd_try   = mean(pd_joint_try, 1);                          % 1 x |remaining|

    % 选 mean_pd 最大的候选（贪婪）
    [best_val, best_local] = max(mean_pd_try);
    best_global = remaining(best_local);

    % 更新状态
    selected(end+1) = best_global;                        %#ok<AGROW>
    remaining(best_local) = [];
    one_minus_cur = one_minus_cur .* (1 - PD_all(:, best_global));

    pd_map_cur     = reshape(1 - one_minus_cur, size(XX));
    mean_pd_hist(step+1)   = best_val;
    coverage_hist(step+1)  = mean(pd_map_cur(:) >= cfg.placement_pd_threshold);
    marginal_gain(step)    = best_val - mean_pd_hist(step);

    if best_val >= cfg.placement_min_coverage
        satisfied = true;
        break;
    end
end

n_nodes = numel(selected);
mean_pd_hist  = mean_pd_hist(1:n_nodes+1);
coverage_hist = coverage_hist(1:n_nodes+1);
marginal_gain = marginal_gain(1:n_nodes);

result = struct( ...
    'sel_idx',       selected, ...
    'sel_coords',    candidates(selected, :), ...
    'mean_pd_hist',  mean_pd_hist, ...
    'coverage_hist', coverage_hist, ...
    'marginal_gain', marginal_gain, ...
    'pd_map_final',  reshape(1 - one_minus_cur, size(XX)), ...
    'xc', xc, 'zc', zc, ...
    'n_nodes',    n_nodes, ...
    'satisfied',  satisfied, ...
    'min_coverage', cfg.placement_min_coverage);

end

% =========================================================================
function [xc, zc] = placement_grid_local(cfg)
if ~isempty(cfg.placement_area_mm)
    ar = cfg.placement_area_mm;
else
    pts = cfg.mic_coords;
    ar  = [min(pts(:,1))-cfg.grid_margin_mm, max(pts(:,1))+cfg.grid_margin_mm, ...
           min(pts(:,3))-cfg.grid_margin_mm, max(pts(:,3))+cfg.grid_margin_mm];
end
res = cfg.placement_area_res_mm;
xc  = ar(1):res:ar(2);
zc  = ar(3):res:ar(4);
end

function Lp = lp_grid_local(mic, A, n, XX, ZZ, cfg)
dx  = XX - mic(1);
dz  = ZZ - mic(3);
dy_d = cfg.src_height_mm - mic(2);
r_d  = max(sqrt(dx.^2 + dy_d.^2 + dz.^2) / 1000.0, 0.05);
Lp_d = A - n .* log10(r_d) - cfg.alpha_atm .* r_d;
if cfg.R_ground > 0
    dy_i = (2*cfg.ground_y_mm - cfg.src_height_mm) - mic(2);
    r_i  = max(sqrt(dx.^2 + dy_i.^2 + dz.^2) / 1000.0, 0.05);
    Lp_i = A - n .* log10(r_i) - cfg.alpha_atm .* r_i;
    Lp   = 10 * log10(10.^(Lp_d/10) + cfg.R_ground.^2 .* 10.^(Lp_i/10));
else
    Lp = Lp_d;
end
end
