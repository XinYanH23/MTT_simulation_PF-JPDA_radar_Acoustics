function seeds = mc_experiment_seeds(n)
%MC_EXPERIMENT_SEEDS  主实验独立、互不重复随机种子。
%
%   seeds = MC_EXPERIMENT_SEEDS()      % 默认 25 个
%   seeds = MC_EXPERIMENT_SEEDS(n)     % n ∈ [20, 30]
%
%   不含调试种子 42。主实验必须用本列表；rng(42) 仅允许 MC_DEBUG=true。

if nargin < 1 || isempty(n)
    n = 25;
end
if n < 1 || n ~= floor(n)
    error('mc_experiment_seeds:Range', 'n_mc 须为正整数。');
end
if n < 20 || n > 30
    warning('mc_experiment_seeds:Range', ...
        '论文主实验建议 20–30 组（默认 25），当前 n=%d。', n);
end

% 1001… 连续整数，与 42 无交集，且彼此不重复
seeds = 1000 + (1:n);
if any(seeds == 42) || numel(unique(seeds)) ~= n
    error('mc_experiment_seeds:Invalid', '种子列表非法。');
end
seeds = seeds(:)';
end
