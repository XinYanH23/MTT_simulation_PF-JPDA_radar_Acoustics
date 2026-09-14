function clusters = jpda_cluster(valid_mat)
% JPDA_CLUSTER  将目标和探测按关联关系划分为不相交的簇。
%
%   clusters = jpda_cluster(valid_mat)
%
%   INPUT
%     valid_mat : T×D logical  验证矩阵（来自 jpda_gating）
%
%   OUTPUT
%     clusters : 1×K cell 数组，clusters{k} 含：
%                  .targets    目标局部索引（对应 tracks 数组下标）
%                  .detections 探测局部索引（对应 Z 列下标）
%
%   算法：
%     在 目标—探测 二部图上做 BFS/DFS 连通分量分解。
%     若目标 t1 和 t2 共享至少一个门内探测，则归入同一簇；
%     反之则独立处理，从而降低联合事件枚举的规模。
%
%   参考：Bar-Shalom (2001) §6.4（Clustering in JPDA）

[T, D] = size(valid_mat);

t_label = zeros(1, T);     % 目标所属簇编号（0 = 未分配）
d_label = zeros(1, D);     % 探测所属簇编号
n_clust = 0;

for t0 = 1 : T
    if t_label(t0) > 0, continue; end          % 已分配
    if ~any(valid_mat(t0, :)), continue; end    % 孤立目标（无门内探测）

    n_clust      = n_clust + 1;
    cid          = n_clust;
    t_label(t0)  = cid;

    % BFS 队列
    q    = t0;
    head = 1;

    while head <= length(q)
        t    = q(head);
        head = head + 1;

        for d = find(valid_mat(t, :))
            if d_label(d) == 0
                d_label(d) = cid;

                % 找所有与探测 d 相关的目标，加入队列
                for t2 = find(valid_mat(:, d))'
                    if t_label(t2) == 0
                        t_label(t2) = cid;
                        q(end + 1)  = t2;   %#ok<AGROW>
                    end
                end
            end
        end
    end
end

% 孤立目标（无任何门内探测）：各自单独成簇
for t0 = 1 : T
    if t_label(t0) == 0
        n_clust     = n_clust + 1;
        t_label(t0) = n_clust;
    end
end

% 构建输出
clusters = cell(1, n_clust);
for k = 1 : n_clust
    clusters{k}.targets    = find(t_label == k);
    clusters{k}.detections = find(d_label == k);
end
end
