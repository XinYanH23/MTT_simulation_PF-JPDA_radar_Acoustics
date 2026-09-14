function p = acoustic_detection_prob(snr_db, cfg)
%ACOUSTIC_DETECTION_PROB 单节点探测概率 P_D(SNR)，sigmoid 形式。
%
%   p = ACOUSTIC_DETECTION_PROB(snr_db, cfg)
%       p = 1 / (1 + exp(-(snr_db - eta) / kappa))
%
%   snr_db 可为标量或数组。参数 eta=det_eta_db, kappa=det_kappa_db。
%   参照 ACOUSTIC_PROBABILISTIC_OBSERVATION_LAYER_RIGOROUS.md §2。

p = 1.0 ./ (1.0 + exp(-(snr_db - cfg.det_eta_db) ./ cfg.det_kappa_db));

end
