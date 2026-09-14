function Theta = jpda_enumerate(T_loc, D_loc, valid_loc, max_events)
% JPDA_ENUMERATE  枚举一个簇内所有可行联合关联事件。
%
%   Theta = jpda_enumerate(T_loc, D_loc, valid_loc, max_events)
%
%   INPUTS
%     T_loc     : 簇内目标数
%     D_loc     : 簇内探测数
%     valid_loc : T_loc×D_loc logical  本簇的局部验证子矩阵
%     max_events: 枚举上限（默认 5000，超出则中止并警告）
%
%   OUTPUT
%     Theta : 1×E cell 数组，每元素为长度 T_loc 的行向量。
%             θ(t) = 0  → 目标 t 本帧未被检测
%             θ(t) = d  → 目标 t 与探测 d 关联（d ∈ 1..D_loc）
%
%   可行性约束（JPDA 标准，Bar-Shalom 2001 §6.2.2）：
%     1. 每个探测最多被一个目标关联
%     2. θ(t) ≠ 0 要求 valid_loc(t, θ(t)) = true
%
%   枚举策略：
%     对目标按顺序递归展开可能的赋值，利用"已占用探测集"剪枝。
%     时间复杂度为 O((D_loc+1)^T_loc) 最坏情况，
%     但在实际稀疏关联下远低于此上界。

if nargin < 4
    max_events = 5000;
end

Theta  = {};
used   = false(1, D_loc);     % 当前路径中已占用的探测

recurse(1, zeros(1, T_loc), used);

    % ── 嵌套递归函数 ─────────────────────────────────────────────────
    function recurse(t, theta, used_flag)
        if length(Theta) >= max_events
            return
        end

        if t > T_loc
            Theta{end + 1} = theta;
            return
        end

        % 选项 A：目标 t 本帧漏检
        theta_new        = theta;
        theta_new(t)     = 0;
        recurse(t + 1, theta_new, used_flag);
        if length(Theta) >= max_events, return; end

        % 选项 B：目标 t 与探测 d 关联
        for d = 1 : D_loc
            if valid_loc(t, d) && ~used_flag(d)
                theta_new        = theta;
                theta_new(t)     = d;
                used_new         = used_flag;
                used_new(d)      = true;
                recurse(t + 1, theta_new, used_new);
                if length(Theta) >= max_events, return; end
            end
        end
    end
% ─────────────────────────────────────────────────────────────────────

if length(Theta) >= max_events
    warning('Tracker:jpda_enum_cap', ...
        'jpda_enumerate: 事件数达到上限 %d，结果为截断枚举。', max_events);
end
end
