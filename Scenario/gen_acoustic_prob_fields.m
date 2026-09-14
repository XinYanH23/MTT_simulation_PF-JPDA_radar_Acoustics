function ac_fields = gen_acoustic_prob_fields(truth, time, acousticPos_m, method)
%GEN_ACOUSTIC_PROB_FIELDS  为所有帧生成声学概率场，供 IMM-PF 主循环使用。
%
%   ac_fields = GEN_ACOUSTIC_PROB_FIELDS(truth, time, acousticPos_m)
%   ac_fields = GEN_ACOUSTIC_PROB_FIELDS(truth, time, acousticPos_m, method)
%
%   INPUTS
%     truth         : 真值目标数组（每元素含 .pos3D N×3，单位 m）
%     time          : N×1 时间向量
%     acousticPos_m : Mx3 麦克风节点坐标 [x,y,z]（单位 m）
%     method        : 概率场方法 'softmax'|'intensity'|'nms_gmm'（默认 'softmax'）
%
%   OUTPUT
%     ac_fields : 1×N cell array，每元素为 build_acoustic_prob_field 的输出 pack
%                 若当帧无目标或所有节点漏检，pack.detected = false
%
%   流程（每帧 k）：
%     1. 计算各目标在各节点的 SPL（球面扩散 + 地面反射）
%     2. 非相干能量叠加（多目标 → 单节点总 SPL）
%     3. 叠加高斯量测噪声
%     4. 调用 build_acoustic_prob_field 生成概率场 pack
%
%   坐标系说明：
%     truth.pos3D 单位 m，acoustic_config_from_scenario 转为 mm
%
%   递归贝叶斯在帧间传递 prev_post 以平滑概率场。

if nargin < 4 || isempty(method); method = 'softmax'; end

N = length(time);
Mn = size(acousticPos_m, 1);

% --- 构建声学配置（使用场景麦克风位置）---
% 估算目标飞行高度（用所有目标第一帧均值作为声源高度参考）
h_vals = [];
for n = 1:length(truth)
    if size(truth(n).pos3D, 1) >= 1
        h_vals(end+1) = truth(n).pos3D(1, 2); %#ok<AGROW>
    end
end
if isempty(h_vals)
    src_h_m = 50.0;
else
    src_h_m = mean(h_vals);
end

ac_cfg = acoustic_config_from_scenario(acousticPos_m, src_h_m);

% --- 预分配输出 ---
ac_fields  = cell(1, N);
prev_posts = repmat({[]}, 1, 1);   % 递归贝叶斯历史（单场景单场，多目标共用同一场）
prev_post  = [];

fprintf('  生成声学概率场（%d 帧）...\n', N);

for k = 1:N
    % --- 收集本帧有效目标位置 ---
    src_list = [];
    for n = 1:length(truth)
        if k > size(truth(n).pos3D, 1); continue; end
        p = truth(n).pos3D(k, :);
        if any(isnan(p)); continue; end
        src_list(end+1, :) = p; %#ok<AGROW>
    end

    if isempty(src_list)
        % 无目标：构造空场（detected=false）
        ac_fields{k} = empty_ac_pack(ac_cfg, method);
        continue
    end

    % --- 各节点：多目标非相干 SPL 叠加 ---
    measured_lp = nan(Mn, 1);
    for m = 1:Mn
        I_total = 0;
        for n_src = 1:size(src_list, 1)
            src_xyz_mm = src_list(n_src, :) * 1000.0;   % m → mm
            Lp_n = single_mic_spl(src_xyz_mm, m, ac_cfg);
            I_total = I_total + 10^(Lp_n / 10);
        end
        Lp_total = 10 * log10(max(I_total, 1e-30));

        % 探测概率（基于总 SNR）
        snr_total = Lp_total - ac_cfg.noise_floor_db(m);
        pd = acoustic_detection_prob(snr_total, ac_cfg);

        if rand() <= pd
            % 加测量噪声
            measured_lp(m) = Lp_total + ac_cfg.sigma0_db * randn();
        end
        % 漏检：保持 NaN
    end

    % --- 构建概率场 ---
    pack = build_acoustic_prob_field(measured_lp, ac_cfg, method, prev_post);

    % 递归贝叶斯：传递后验
    prev_post = pack.post;

    ac_fields{k} = pack;

    if mod(k, 100) == 0
        fprintf('    帧 %d/%d\n', k, N);
    end
end

fprintf('  声学概率场生成完毕。\n');
end

% =========================================================================
function pack = empty_ac_pack(ac_cfg, method)
%EMPTY_AC_PACK  无目标帧的空 pack（detected=false，概率场为均匀分布）
if nargin < 2; method = 'softmax'; end

measured_lp = nan(size(ac_cfg.mic_coords, 1), 1);
pack = build_acoustic_prob_field(measured_lp, ac_cfg, method, []);
pack.detected = false;
end
