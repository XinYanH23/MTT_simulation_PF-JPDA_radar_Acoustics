function [pos_pool, extra] = expand_acoustic_catalog(pos_base, buildings, n_max, seed)
%EXPAND_ACOUSTIC_CATALOG  在已有声学站外按最远点补点，得到嵌套站点池。
%
%   pos_pool 前 size(pos_base,1) 行保持原场景 33 站不动；
%   其后为补点。n<=33 仍从原 33 站做 FPS；n>33 为「全部原站 + 补点前缀」。

if nargin < 4 || isempty(seed); seed = 2024; end
M = size(pos_base, 1);
n_max = max(round(n_max), M);
if n_max <= M
    pos_pool = pos_base;
    extra = zeros(0, 3);
    return
end

need = n_max - M;
z0 = median(pos_base(:, 3));
cands = zeros(0, 3);

this_dir = fileparts(mfilename('fullpath'));
v2_dir = fullfile(this_dir, '..', '..', 'v2.0');
if exist(v2_dir, 'dir')
    addpath(v2_dir);
end
if nargin >= 2 && ~isempty(buildings) && exist('generateSensorCandidates', 'file')
    try
        cands = [cands; generateSensorCandidates(buildings, 1.0)]; %#ok<AGROW>
    catch
    end
end

[xg, yg] = meshgrid(25:50:975, 25:50:975);
grid_xyz = [xg(:), yg(:), z0 * ones(numel(xg), 1)];
cands = [cands; grid_xyz];

% 丢掉已有站 35 m 内的候选，避免叠在原站上
keep = true(size(cands, 1), 1);
for i = 1:size(cands, 1)
    dmin = min(sqrt(sum((pos_base(:, 1:2) - cands(i, 1:2)).^2, 2)));
    if dmin < 35
        keep(i) = false;
    end
end
cands = unique(cands(keep, :), 'rows');
if size(cands, 1) < need
    error('expand_acoustic_catalog:Short', ...
        '候选不足：需要 %d 个补点，只剩 %d。', need, size(cands, 1));
end

rng(seed);
extra = zeros(need, 3);
selected_xy = pos_base(:, 1:2);
dmin = min_dist2(cands(:, 1:2), selected_xy);
for k = 1:need
    [~, ik] = max(dmin);
    extra(k, :) = cands(ik, :);
    selected_xy = [selected_xy; cands(ik, 1:2)]; %#ok<AGROW>
    dmin = min(dmin, sum((cands(:, 1:2) - cands(ik, 1:2)).^2, 2));
    dmin(ik) = -1;
end

pos_pool = [pos_base; extra];
end

function d2 = min_dist2(xy, ref)
d2 = inf(size(xy, 1), 1);
for j = 1:size(ref, 1)
    d2 = min(d2, sum((xy - ref(j, :)).^2, 2));
end
end
