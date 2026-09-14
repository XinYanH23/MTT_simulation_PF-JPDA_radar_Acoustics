function diag = diagnose_track_redundancy(result_mat, assoc_radius)
%DIAGNOSE_TRACK_REDUNDANCY  诊断"5 个真目标 → 56 条确认轨迹"的成因。
%
%   diag = DIAGNOSE_TRACK_REDUNDANCY()              % 读默认结果文件
%   diag = DIAGNOSE_TRACK_REDUNDANCY(result_mat)    % 指定 .mat
%   diag = DIAGNOSE_TRACK_REDUNDANCY(result_mat, R) % 自定义贴合半径(m)
%
%   区分三种成因：
%     冗余(redundancy)   ：同一真目标同一帧被 ≥2 条确认轨迹同时覆盖
%     碎裂(fragmentation)：同一真目标整段被多条确认轨迹"接力"覆盖
%     杂波(clutter)      ：确认轨迹不贴任何真目标
%
%   不修改 Tracker，任何阶段只读 Step3 结果。

if nargin < 1 || isempty(result_mat)
    here       = fileparts(mfilename('fullpath'));
    result_mat = fullfile(here, '..', 'Step3_IMM_JPDA_PF_Result.mat');
end
if nargin < 2 || isempty(assoc_radius)
    assoc_radius = 30;   % m，判定"该确认轨迹贴在某真目标上"的距离阈值
end

S = load(result_mat, 'tracks', 'truth', 'time', 'cfg');
tracks = S.tracks; truth = S.truth; time = S.time; cfg = S.cfg;
N = numel(time); Nt = numel(truth);

% ── 取"曾确认"轨迹（与 eval_metrics 口径一致）──────────────────────────
conf = false(1, numel(tracks));
for i = 1:numel(tracks)
    conf(i) = strcmp(tracks(i).status,'confirmed') || tracks(i).total_hits >= cfg.mn_M;
end
conf_idx = find(conf);
nc = numel(conf_idx);

% ── 逐帧：每个真目标被多少条确认轨迹覆盖 ───────────────────────────────
cover_cnt   = zeros(Nt, N);          % 真目标 j 在帧 k 被几条确认轨迹覆盖
est_per_fr  = zeros(1, N);           % 帧 k 的确认估计数
truth_per_fr= zeros(1, N);           % 帧 k 的活跃真值数
track_hits_truth = false(nc, Nt);    % 确认轨迹 ci 是否曾贴过真目标 j

for k = 1:N
    % 活跃真值
    gt = nan(2, Nt);
    for j = 1:Nt
        raw = truth(j).pos3D;
        if k <= size(raw,1) && ~any(isnan(raw(k,1:2)))
            gt(:,j) = raw(k,1:2)';
        end
    end
    truth_per_fr(k) = sum(~isnan(gt(1,:)));

    % 本帧确认估计
    for ci = 1:nc
        i = conf_idx(ci);
        th = tracks(i).time_hist;
        [dmin, loc] = min(abs(th - k));
        if isempty(loc) || dmin > 1, continue; end
        est = tracks(i).state_hist(1:2, loc);
        est_per_fr(k) = est_per_fr(k) + 1;
        for j = 1:Nt
            if isnan(gt(1,j)), continue; end
            if norm(est - gt(:,j)) <= assoc_radius
                cover_cnt(j,k) = cover_cnt(j,k) + 1;
                track_hits_truth(ci,j) = true;
            end
        end
    end
end

% ── 指标汇总 ───────────────────────────────────────────────────────────
active = truth_per_fr > 0;
% 1) 冗余因子：活跃帧里，平均每个"被覆盖"的真目标摊到几条确认轨迹
cov_active = cover_cnt(:, active);
covered = cov_active(cov_active > 0);
diag.redundancy_factor   = mean(covered);                 % 期望≈1，越大越冗余
diag.max_concurrent_per_truth = max(cover_cnt(:));
% 2) 估计/真值 基数比（OSPA 基数误差来源）
diag.mean_est_per_frame   = mean(est_per_fr(active));
diag.mean_truth_per_frame = mean(truth_per_fr(active));
diag.cardinality_ratio    = diag.mean_est_per_frame / max(diag.mean_truth_per_frame,eps);
% 3) 碎裂：每个真目标整段被多少条不同确认轨迹覆盖过
frag = sum(track_hits_truth, 1);                          % 1×Nt
diag.fragmentation_per_truth = frag;
diag.mean_fragmentation = mean(frag);                     % 越大越碎裂/冗余
% 4) 杂波：不贴任何真目标的确认轨迹数
clutter = sum(~any(track_hits_truth, 2));
diag.clutter_tracks = clutter;
diag.clutter_ratio  = clutter / max(nc,1);
diag.n_confirmed    = nc;
diag.n_truth        = Nt;

% ── 打印 ───────────────────────────────────────────────────────────────
fprintf('\n===== 轨迹冗余诊断（半径 %.0f m）=====\n', assoc_radius);
fprintf('  曾确认轨迹数 / 真目标数      : %d / %d\n', nc, Nt);
fprintf('  平均估计数/帧 vs 真值数/帧   : %.1f vs %.1f  (基数比 %.1fx)\n', ...
        diag.mean_est_per_frame, diag.mean_truth_per_frame, diag.cardinality_ratio);
fprintf('  冗余因子(每被覆盖目标的轨迹数): %.2f  (>1 即并发冗余)\n', diag.redundancy_factor);
fprintf('  单目标最大并发确认轨迹数      : %d\n', diag.max_concurrent_per_truth);
fprintf('  碎裂因子(每真目标涉及轨迹数)  : %.2f  [', diag.mean_fragmentation);
fprintf('%d ', frag); fprintf(']\n');
fprintf('  纯杂波确认轨迹(不贴任何真值)  : %d / %d  (%.0f%%)\n', ...
        clutter, nc, 100*diag.clutter_ratio);
fprintf('---------------------------------------------\n');
if diag.redundancy_factor > 1.5
    fprintf('  ⇒ 主因：并发冗余（同目标多条确认轨迹并存）。\n');
elseif diag.mean_fragmentation > 2 && diag.redundancy_factor <= 1.5
    fprintf('  ⇒ 主因：轨迹碎裂（同目标接力重建）。\n');
end
if diag.clutter_ratio > 0.2
    fprintf('  ⇒ 兼有：杂波假轨迹占比偏高。\n');
end
fprintf('=============================================\n\n');
end
