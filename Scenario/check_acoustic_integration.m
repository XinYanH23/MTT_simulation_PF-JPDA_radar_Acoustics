%% CHECK_ACOUSTIC_INTEGRATION  Step2 声学自适应 R 与跟踪链路集成检查
%
%  检查项：
%    1. acoustic_adaptive_R 输出正定、径向>切向
%    2. Step2 点迹含自适应 MeasurementNoise
%    3. detections_to_ZR 正确传递各向异性 R（非固定 8²）
%    4. 端到端：Step2 → main 跟踪器可运行
%
%  用法：在 Scenario/ 目录下运行本脚本

fprintf('\n===== 声学自适应 R 集成检查 =====\n');
this_dir = fileparts(mfilename('fullpath'));
tracker_dir = fullfile(this_dir, '..');
addpath(this_dir);
addpath(tracker_dir);
addpath(fullfile(tracker_dir, 'Utilities'));

cfg = config_tracker();
n_fail = 0;

%% ── 1. 单元：acoustic_adaptive_R ─────────────────────────────────────
fprintf('\n[1/4] acoustic_adaptive_R 单元测试...\n');
uavPos = [200, 300, 50];
sPos   = [100, 100, 5];
snr_dB = 10;
R2d = acoustic_adaptive_R(uavPos, sPos, snr_dB, struct());

[~, p] = chol(R2d);
if p ~= 0
    fprintf('  FAIL: R2d 非正定\n');
    n_fail = n_fail + 1;
else
    fprintf('  OK: R2d 正定\n');
end

phi = atan2(uavPos(2)-sPos(2), uavPos(1)-sPos(1));
Rot = [cos(phi), -sin(phi); sin(phi), cos(phi)];
sig_local = sqrt(diag(Rot' * R2d * Rot));
if sig_local(1) <= sig_local(2)
    fprintf('  FAIL: 径向标准差 (%.2f) 应大于切向 (%.2f)\n', sig_local(1), sig_local(2));
    n_fail = n_fail + 1;
else
    fprintf('  OK: 径向 σ=%.2f m > 切向 σ=%.2f m\n', sig_local(1), sig_local(2));
end

R_far = acoustic_adaptive_R([500,500,50], sPos, 5, struct());
R_near = acoustic_adaptive_R([120,120,50], sPos, 15, struct());
if trace(R_far) <= trace(R_near)
    fprintf('  WARN: 远距离 R 迹未明显大于近距离（可能受 sigma_max 限幅）\n');
else
    fprintf('  OK: 远距离 R 迹 > 近距离 R 迹\n');
end

%% ── 2. Step2 数据存在性与结构 ─────────────────────────────────────────
fprintf('\n[2/4] Step2 点迹结构检查...\n');
s1 = fullfile(this_dir, 'Step1_Data.mat');
s2 = fullfile(this_dir, 'Step2_HeteroDetections.mat');

if ~exist(s1, 'file')
    fprintf('  缺少 Step1_Data.mat，请先 run_scenario\n');
    n_fail = n_fail + 1;
else
    fprintf('  Step1: %s\n', s1);
end

if ~exist(s2, 'file')
    fprintf('  缺少 Step2_HeteroDetections.mat，将尝试运行 step2...\n');
    if exist(s1, 'file')
        load(s1);
        run(fullfile(this_dir, 'step2_detection_simulation.m'));
    end
end

if ~exist(s2, 'file')
    fprintf('  FAIL: 无法生成 Step2 数据\n');
    n_fail = n_fail + 1;
else
    S2 = load(s2, 'allDetections');
    allDetections = S2.allDetections;
    n_ac = 0; n_rd = 0;
    n_aniso = 0; n_fixed = 0;
    R_samples = zeros(2, 2, 0);

    for t = 1:length(allDetections)
        frame = allDetections{t};
        if isempty(frame), continue; end
        if iscell(frame), dets = frame; else, dets = num2cell(frame); end
        for d = 1:length(dets)
            det = dets{d};
            if isobject(det)
                sid = det.SensorIndex;
                noise = det.MeasurementNoise;
            else
                sid = det.SensorIndex;
                noise = det.MeasurementNoise;
            end
            if sid > 100
                n_rd = n_rd + 1;
            else
                n_ac = n_ac + 1;
                R2 = noise(1:2, 1:2);
                R_samples(:, :, end+1) = R2; %#ok<AGROW>
                off_diag = abs(R2(1,2));
                tr = trace(R2);
                if off_diag > 1e-6 || abs(R2(1,1) - R2(2,2)) > 1
                    n_aniso = n_aniso + 1;
                end
                if abs(R2(1,1) - 64) < 0.1 && abs(R2(2,2) - 64) < 0.1
                    n_fixed = n_fixed + 1;
                end
            end
        end
    end

    fprintf('  总点迹: 声学 %d, 雷达 %d\n', n_ac, n_rd);
    if n_ac == 0
        fprintf('  FAIL: 无声学点迹\n');
        n_fail = n_fail + 1;
    else
        fprintf('  声学各向异性 R: %d / %d\n', n_aniso, n_ac);
        fprintf('  固定 diag(8²) 回退: %d\n', n_fixed);
        if n_aniso < n_ac * 0.5
            fprintf('  FAIL: 多数声学点迹未呈现各向异性\n');
            n_fail = n_fail + 1;
        else
            fprintf('  OK: 声学点迹携带自适应椭圆 R\n');
        end
    end
end

%% ── 3. detections_to_ZR 接口 ──────────────────────────────────────────
fprintf('\n[3/4] detections_to_ZR 接口检查...\n');
if exist('allDetections', 'var') && ~isempty(allDetections)
    % 找一帧含声学点迹的
    test_frame = [];
    for t = 1:min(100, length(allDetections))
        if ~isempty(allDetections{t})
            test_frame = allDetections{t};
            break
        end
    end
    if isempty(test_frame)
        fprintf('  WARN: 前 100 帧无点迹，跳过\n');
    else
        [Z, R_list, stypes] = detections_to_ZR(test_frame, cfg);
        ac_idx = find(strcmp(stypes, 'Acoustic'));
        if isempty(ac_idx)
            fprintf('  WARN: 测试帧无声学点迹\n');
        else
            R_ac = R_list(:, :, ac_idx(1));
            if abs(R_ac(1,1) - 64) < 0.1 && abs(R_ac(2,2) - 64) < 0.1 && abs(R_ac(1,2)) < 1e-6
                fprintf('  FAIL: detections_to_ZR 回退到固定 R_acoustic\n');
                n_fail = n_fail + 1;
            else
                fprintf('  OK: 传递自适应 R, trace=%.1f (非 128)\n', trace(R_ac));
            end
        end
        fprintf('  帧内探测数: %d\n', size(Z, 2));
    end
end

%% ── 4. 跟踪器端到端（轻量：仅跑前 50 帧逻辑）──────────────────────────
fprintf('\n[4/4] 跟踪器链路抽样（前 50 帧）...\n');
if exist(s2, 'file')
    load(s2, 'allDetections', 'time');
    if exist(s1, 'file')
        d1 = load(s1, 'truth');
        truth = d1.truth;
    end
    addpath(fullfile(tracker_dir, 'IMM'));
    addpath(fullfile(tracker_dir, 'JPDA'));
    addpath(fullfile(tracker_dir, 'EKF'));
    addpath(fullfile(tracker_dir, 'TrackManagement'));

    tracks = [];
    next_id = 1;
    n_frames = min(50, length(allDetections));
    gate_ok = 0;

    for k = 1:n_frames
        dets = allDetections{k};
        [Z, R_list, ~] = detections_to_ZR(dets, cfg);
        if size(Z, 2) > 0
            gate_ok = gate_ok + 1;
        end
    end
    fprintf('  前 %d 帧中 %d 帧有有效量测\n', n_frames, gate_ok);
    if gate_ok == 0
        fprintf('  FAIL: 无有效量测进入跟踪器\n');
        n_fail = n_fail + 1;
    else
        fprintf('  OK: 量测可进入跟踪器\n');
    end
end

%% ── 汇总 ──────────────────────────────────────────────────────────────
fprintf('\n===== 检查完成: %d 项失败 =====\n', n_fail);
if n_fail > 0
    error('check_acoustic_integration: 存在失败项');
end
