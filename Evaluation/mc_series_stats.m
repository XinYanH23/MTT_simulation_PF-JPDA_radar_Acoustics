function s = mc_series_stats(X)
%MC_SERIES_STATS  对 n_mc×N 矩阵逐列做均值 / 标准差 / 95% CI。
[n, N] = size(X);
s.n     = n;
s.mean  = mean(X, 1, 'omitnan');
if n <= 1
    s.std = zeros(1, N);
else
    s.std = std(X, 0, 1, 'omitnan');
end
s.sem   = s.std / sqrt(max(n, 1));
tcrit   = mc_tinv_975(max(n - 1, 1));
s.ci_lo = s.mean - tcrit * s.sem;
s.ci_hi = s.mean + tcrit * s.sem;
end
