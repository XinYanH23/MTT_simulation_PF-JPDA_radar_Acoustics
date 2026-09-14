function [seeds, n_mc, is_debug] = mc_resolve_seeds(MC_DEBUG, n_mc)
%MC_RESOLVE_SEEDS  调试用 rng(42)；主实验用 25 组独立种子。
if nargin < 1 || isempty(MC_DEBUG)
    MC_DEBUG = false;
end
if nargin < 2 || isempty(n_mc)
    n_mc = 25;
end
if MC_DEBUG
    seeds = 42;
    n_mc = 1;
    is_debug = true;
    fprintf(['【调试】rng(42) 单次运行（仅用于排错；' ...
             '主实验须 MC_DEBUG=false，25 组独立种子）\n']);
else
    seeds = mc_experiment_seeds(n_mc);
    is_debug = false;
    fprintf('【主实验】%d 组独立种子（不含 42）：%s\n', ...
        numel(seeds), mat2str(seeds));
end
end
