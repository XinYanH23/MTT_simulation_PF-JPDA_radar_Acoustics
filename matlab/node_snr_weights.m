function [snr_db, w, p_detect_m, p_detect_joint] = node_snr_weights(measured_lp, cfg)
%NODE_SNR_WEIGHTS 由各节点实测 SPL 计算 SNR、节点可靠性权重与探测概率。
%
%   [snr_db, w, p_detect_m, p_detect_joint] = NODE_SNR_WEIGHTS(measured_lp, cfg)
%
%   输入：
%     measured_lp : Mx1 各节点实测声压级 (dB)，缺测节点可填 NaN
%   输出：
%     snr_db         : Mx1 各节点 SNR = measured_lp - noise_floor
%     w              : Mx1 节点可靠性权重，softmax(snr/weight_tau)，缺测置0，归一化使 sum(w)=有效节点数
%     p_detect_m     : Mx1 各节点探测概率
%     p_detect_joint : 标量，至少一节点探测到 = 1 - prod(1 - p_detect_m)
%
%   这是"多频带提取"在无音频仿真下的等效：以节点级 SNR 作为可靠性，
%   高 SNR 节点在跨节点对数域融合中获得更大权重（参照 4SPL/band_weighting.py 思想）。

measured_lp = measured_lp(:);
M = numel(measured_lp);
nf = cfg.noise_floor_db(:);

snr_db = measured_lp - nf;
valid = ~isnan(snr_db);

% 探测概率
p_detect_m = zeros(M, 1);
p_detect_m(valid) = acoustic_detection_prob(snr_db(valid), cfg);
p_detect_joint = 1 - prod(1 - p_detect_m(valid));
if isempty(p_detect_m(valid)); p_detect_joint = 0; end

% 节点可靠性权重：softmax(snr / tau)
w = zeros(M, 1);
if any(valid)
    s = snr_db(valid);
    s = s - max(s);                         % 数值稳定
    e = exp(s / cfg.weight_tau);
    wv = e / sum(e);
    wv = wv * sum(valid);                   % 归一化使平均权重=1（保持 chi2 量纲）
    w(valid) = wv;
end

end
