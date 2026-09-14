function str = mc_pm_tex(s, kind)
%MC_PM_TEX  均值 ± 标准差 的 LaTeX 单元格。
if nargin < 2 || isempty(kind)
    kind = 'm';
end
switch lower(kind)
    case 'pct'
        str = sprintf('%.1f $\\pm$ %.1f\\%%', 100*s.mean, 100*s.std);
    case 'f1'
        str = sprintf('%.1f $\\pm$ %.1f', s.mean, s.std);
    otherwise
        str = sprintf('%.2f $\\pm$ %.2f', s.mean, s.std);
end
end
