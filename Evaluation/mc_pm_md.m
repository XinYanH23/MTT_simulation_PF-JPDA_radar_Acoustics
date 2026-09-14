function str = mc_pm_md(s, kind)
%MC_PM_MD  均值 ± 标准差（Markdown / 控制台）。
if nargin < 2 || isempty(kind)
    kind = 'm';
end
switch lower(kind)
    case 'pct'
        str = sprintf('%.1f ± %.1f%%', 100*s.mean, 100*s.std);
    case 'f1'
        str = sprintf('%.1f ± %.1f', s.mean, s.std);
    otherwise
        str = sprintf('%.2f ± %.2f', s.mean, s.std);
end
end
