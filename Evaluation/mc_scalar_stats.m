function s = mc_scalar_stats(x)
%MC_SCALAR_STATS  样本均值、标准差、95% t 置信区间。
x = x(:);
x = x(isfinite(x));
n = numel(x);
s = struct('n', n, 'mean', NaN, 'std', NaN, 'sem', NaN, ...
           'ci_lo', NaN, 'ci_hi', NaN);
if n == 0
    return
end
s.mean = mean(x);
if n == 1
    s.std = 0;
    s.sem = 0;
    s.ci_lo = s.mean;
    s.ci_hi = s.mean;
    return
end
s.std   = std(x, 0);
s.sem   = s.std / sqrt(n);
tcrit   = mc_tinv_975(n - 1);
s.ci_lo = s.mean - tcrit * s.sem;
s.ci_hi = s.mean + tcrit * s.sem;
end
