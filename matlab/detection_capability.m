function out = detection_capability(mic_idx, cfg, r_axis_m)
%DETECTION_CAPABILITY 单节点声学探测能力随距离的变化曲线。
%
%   out = DETECTION_CAPABILITY(mic_idx, cfg)
%   out = DETECTION_CAPABILITY(mic_idx, cfg, r_axis_m)
%
%   沿径向（声源正上方等效）扫描距离 r，计算：
%     Lp(r)   单点 SPL（含地面反射）
%     SNR(r)  = Lp(r) - noise_floor
%     P_D(r)  探测概率
%   并求出 P_D=0.5 对应的最大可探测距离 r50。
%
%   输出 struct：r_m, Lp_db, snr_db, p_detect, r50_m
%
%   注：距离扫描以"声源在麦克风径向外、保持 cfg.src_height_mm 高度"近似，
%       用于刻画该节点的探测能力包络，非某条具体轨迹。

if nargin < 3 || isempty(r_axis_m)
    r_axis_m = linspace(0.5, 60, 200);
end

mic = cfg.mic_coords(mic_idx, :);
nf  = cfg.noise_floor_db(mic_idx);

% 沿 +X 方向在固定高度放置声源，水平距离为 r 的水平分量
h = cfg.src_height_mm - mic(2);   % 垂直高度差 (mm)
Lp_db = zeros(size(r_axis_m));
for k = 1:numel(r_axis_m)
    r_mm = r_axis_m(k) * 1000.0;
    % 让斜距等于 r：水平分量 = sqrt(r^2 - h^2) 若 r>h，否则贴近正上方
    horiz = sqrt(max(r_mm.^2 - h.^2, 0));
    src = [mic(1) + horiz, cfg.src_height_mm, mic(3)];
    Lp_db(k) = single_mic_spl(src, mic_idx, cfg);
end

snr_db = Lp_db - nf;
p_detect = acoustic_detection_prob(snr_db, cfg);

% r50：P_D 首次降到 0.5 以下的距离（探测能力边界）
r50_m = NaN;
below = find(p_detect < 0.5, 1, 'first');
if ~isempty(below) && below > 1
    % 线性插值更精确
    x1 = r_axis_m(below-1); x2 = r_axis_m(below);
    y1 = p_detect(below-1); y2 = p_detect(below);
    r50_m = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
elseif ~isempty(below) && below == 1
    r50_m = r_axis_m(1);
end

out = struct('r_m', r_axis_m, 'Lp_db', Lp_db, 'snr_db', snr_db, ...
             'p_detect', p_detect, 'r50_m', r50_m, 'mic_idx', mic_idx);

end
