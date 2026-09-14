function pf = pf_systematic_resample(pf)
%PF_SYSTEMATIC_RESAMPLE  系统重采样，O(N)，重采样后权重均匀化为 1/N。
%
%   参考：Gustafsson (2010), "Particle Filter Theory and Practice with
%         Positioning Applications", Algorithm 2.

N = pf.N;
w = pf.weights(:)';
C = cumsum(w);  C(end) = 1.0;            % 防浮点误差

u1  = rand() / N;
U   = u1 + (0:N-1) / N;                 % 1 x N 等间距格点

idx = zeros(1, N);
j   = 1;
for m = 1:N
    while j < N && C(j) < U(m);  j = j+1;  end
    idx(m) = j;
end

pf.particles = pf.particles(:, idx);
pf.weights   = ones(1, N) / N;
end
