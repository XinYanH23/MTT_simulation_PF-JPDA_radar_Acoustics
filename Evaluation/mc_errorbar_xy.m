function h = mc_errorbar_xy(x, stats_list, color, marker)
%MC_ERRORBAR_XY  横轴离散点上的均值 ± 95% CI（扫描图）。
if nargin < 4 || isempty(marker)
    marker = 'o';
end
x = x(:)';
n = numel(stats_list);
mu = zeros(1, n); lo = zeros(1, n); hi = zeros(1, n);
for i = 1:n
    mu(i) = stats_list{i}.mean;
    lo(i) = stats_list{i}.ci_lo;
    hi(i) = stats_list{i}.ci_hi;
end
neg = max(mu - lo, 0);
pos = max(hi - mu, 0);
h = errorbar(x, mu, neg, pos, ['-' marker], ...
    'Color', color, 'LineWidth', 1.5, 'MarkerFaceColor', color, 'CapSize', 5);
end
