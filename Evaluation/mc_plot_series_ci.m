function h = mc_plot_series_ci(t, ser, color, stride)
%MC_PLOT_SERIES_CI  均值曲线 + 95% CI 误差棒（稀疏取样以免过密）。
if nargin < 4 || isempty(stride)
    stride = 25;
end
t = t(:)';
h = plot(t, ser.mean, '-', 'Color', color, 'LineWidth', 1.3);
hold on;
idx = unique([1:stride:numel(t), numel(t)]);
neg = ser.mean(idx) - ser.ci_lo(idx);
pos = ser.ci_hi(idx) - ser.mean(idx);
neg = max(neg, 0);
pos = max(pos, 0);
errorbar(t(idx), ser.mean(idx), neg, pos, ...
    'Color', color, 'LineStyle', 'none', 'LineWidth', 0.9, 'CapSize', 3);
end
