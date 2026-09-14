function [pos_out, idx] = subsample_acoustic_stations(pos_in, n_keep, seed)
%SUBSAMPLE_ACOUSTIC_STATIONS  贪心最远点采样，得到稀疏声学站子集。
%
%   [pos_out, idx] = subsample_acoustic_stations(pos_in, n_keep)
%   [pos_out, idx] = subsample_acoustic_stations(pos_in, n_keep, seed)
%
%   pos_in : Mx3 站点坐标；n_keep <= M 时做 FPS，否则原样返回。

if nargin < 3 || isempty(seed); seed = 42; end
M = size(pos_in, 1);
n_keep = min(max(round(n_keep), 1), M);
if n_keep >= M
    pos_out = pos_in;
    idx = (1:M)';
    return
end

rng(seed);
xy = pos_in(:, 1:2);
% 起点：最靠近场景中心的站，保证可复现且覆盖中心
c = mean(xy, 1);
[~, i0] = min(sum((xy - c).^2, 2));
idx = zeros(n_keep, 1);
idx(1) = i0;
dmin = sum((xy - xy(i0,:)).^2, 2);
for k = 2:n_keep
    [~, ik] = max(dmin);
    idx(k) = ik;
    dmin = min(dmin, sum((xy - xy(ik,:)).^2, 2));
end
idx = sort(idx);
pos_out = pos_in(idx, :);
end
