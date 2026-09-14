function tracks = track_merge(tracks, cfg)
%TRACK_MERGE  并发冗余轨迹去重（每帧在 track_manage 之后调用一次）。
%
%   tracks = TRACK_MERGE(tracks, cfg)
%
%   现象：雷达漏关联 + 初始门宽 + birth_vel=0，导致同一真目标被多条轨迹
%   同时覆盖（并发冗余），OSPA 基数误差爆表。本函数对活跃轨迹两两比较，若
%       位置距离 <= cfg.merge_radius  且  速度差满足门限
%   则判为同一目标的重复轨迹，保留质量更高者，其余标记为 deleted。
%
%   速度门：
%     - 两条航迹年龄均 >= merge_young_age：要求 ||dv|| <= merge_vel
%     - 任一条更年轻：仅用位置门（新生航迹速度常为 0，速度门会失效）
%   可选合并 tentative（cfg.merge_tentative=true），避免确认前堆积冗余。

if isempty(tracks) || ~isstruct(tracks)
    return
end
if isfield(cfg, 'track_merge_enable') && ~cfg.track_merge_enable
    return
end

radius = cfg.merge_radius;
vgate  = cfg.merge_vel;
young  = 15;
if isfield(cfg, 'merge_young_age') && ~isempty(cfg.merge_young_age)
    young = cfg.merge_young_age;
end
incl_tent = isfield(cfg, 'merge_tentative') && cfg.merge_tentative;

% 活跃 confirmed（可选 + tentative）
keep = false(1, numel(tracks));
for i = 1:numel(tracks)
    st = tracks(i).status;
    if strcmp(st, 'confirmed')
        keep(i) = true;
    elseif incl_tent && strcmp(st, 'tentative')
        keep(i) = true;
    end
end
idx = find(keep);
n   = numel(idx);
if n < 2
    return
end

for a = 1 : n
    ia = idx(a);
    if strcmp(tracks(ia).status, 'deleted'), continue; end
    for b = a+1 : n
        ib = idx(b);
        if strcmp(tracks(ib).status, 'deleted'), continue; end

        dp = norm(tracks(ia).x(1:2) - tracks(ib).x(1:2));
        if dp > radius
            continue
        end

        % 年轻航迹：速度估计不可靠，跳过速度门
        age_a = tracks(ia).age;
        age_b = tracks(ib).age;
        if age_a >= young && age_b >= young
            dv = norm(tracks(ia).x(3:4) - tracks(ib).x(3:4));
            if dv > vgate
                continue
            end
        end

        % 保留质量更高者；tentative 永远让位于 confirmed
        if prefer(tracks(ib), tracks(ia))
            tracks(ia).status = 'deleted';
            ia = ib;
        else
            tracks(ib).status = 'deleted';
        end
    end
end
end

% -------------------------------------------------------------------------
function tf = prefer(a, b)
%PREFER  true 表示 a 优于 b（应保留 a）。
%  优先级：confirmed > tentative；命中数；更年长；更小 ID。
conf_a = strcmp(a.status, 'confirmed');
conf_b = strcmp(b.status, 'confirmed');
if conf_a ~= conf_b
    tf = conf_a;
    return
end
qa = a.total_hits * 1e6 - a.birth_frame * 1e1 - double(a.id) * 1e-3;
qb = b.total_hits * 1e6 - b.birth_frame * 1e1 - double(b.id) * 1e-3;
tf = qa > qb;
end
