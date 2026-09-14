function mc_bar_with_ci(stats_list, colors, tick_labels)
%MC_BAR_WITH_CI  柱状均值 + 95% CI 误差棒。
n = numel(stats_list);
mu = zeros(n, 1); lo = zeros(n, 1); hi = zeros(n, 1);
for i = 1:n
    mu(i) = stats_list{i}.mean;
    lo(i) = stats_list{i}.ci_lo;
    hi(i) = stats_list{i}.ci_hi;
end
b = bar(mu, 0.58);
b.FaceColor = 'flat';
for i = 1:n
    b.CData(i, :) = colors(i, :);
end
hold on;
neg = max(mu - lo, 0);
pos = max(hi - mu, 0);
errorbar(1:n, mu, neg, pos, 'k.', 'LineWidth', 1.0, 'CapSize', 8);
if nargin >= 3 && ~isempty(tick_labels)
    set(gca, 'XTick', 1:n, 'XTickLabel', tick_labels);
end
end
