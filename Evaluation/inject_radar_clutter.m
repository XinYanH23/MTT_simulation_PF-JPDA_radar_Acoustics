function [allDets_out, stats] = inject_radar_clutter(allDets_in, opts)
%INJECT_RADAR_CLUTTER  向每帧雷达量测集注入均匀泊松杂波点迹。
%
%   opts.mu_fa_extra   : 每帧额外虚警均值（泊松）
%   opts.map_xlim      : [xmin xmax]，默认 [0 1000]
%   opts.map_ylim      : [ymin ymax]，默认 [0 1000]
%   opts.R_radar       : 2×2 量测噪声，默认 diag([1 1])
%   opts.seed          : RNG 种子（可选；不设则不改当前 RNG）
%   opts.z_height      : 杂波高度，默认 50 m（仅填充 Measurement 第 3 维）

if nargin < 2 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'mu_fa_extra'); opts.mu_fa_extra = 0; end
if ~isfield(opts, 'map_xlim');    opts.map_xlim = [0, 1000]; end
if ~isfield(opts, 'map_ylim');    opts.map_ylim = [0, 1000]; end
if ~isfield(opts, 'R_radar');     opts.R_radar = diag([1, 1]); end
if ~isfield(opts, 'z_height');    opts.z_height = 50; end

if isfield(opts, 'seed') && ~isempty(opts.seed)
    rng(opts.seed);
end

N = numel(allDets_in);
allDets_out = allDets_in;
n_injected = zeros(N, 1);
R3 = eye(3);
R3(1:2,1:2) = opts.R_radar(1:2,1:2);

xl = opts.map_xlim; yl = opts.map_ylim;
mu = max(opts.mu_fa_extra, 0);

for k = 1:N
    n_add = poissrnd(mu);
    n_injected(k) = n_add;
    if n_add == 0, continue; end

    frame = allDets_out{k};
    if isempty(frame)
        frame = {};
    elseif ~iscell(frame)
        frame = num2cell(frame);
    end

    for c = 1:n_add
        xy = [xl(1) + (xl(2)-xl(1))*rand; ...
              yl(1) + (yl(2)-yl(1))*rand];
        det = struct();
        det.Measurement = [xy; opts.z_height];
        det.MeasurementNoise = R3;
        det.ObjectAttributes = struct( ...
            'ClassConfidence', 0.1, ...
            'SensorType', 'Radar', ...
            'IsInjectedClutter', true);
        frame{end+1} = det; %#ok<AGROW>
    end
    allDets_out{k} = frame;
end

stats = struct();
stats.mu_fa_extra = mu;
stats.n_injected_total = sum(n_injected);
stats.n_injected_mean = mean(n_injected);
stats.n_injected_per_frame = n_injected;
stats.map_area = (xl(2)-xl(1)) * (yl(2)-yl(1));
stats.lambda_spatial = stats.n_injected_mean / max(stats.map_area, eps);
end
