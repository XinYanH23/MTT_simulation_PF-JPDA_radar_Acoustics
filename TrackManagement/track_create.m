function track = track_create(z, R, track_id, frame_idx, cfg)
% TRACK_CREATE  从一个未关联探测创建新的 Tentative 轨迹。
%
%   track = track_create(z, R, track_id, frame_idx, cfg)
%
%   INPUTS
%     z         : 2×1  位置量测 [px py]'
%     R         : 2×2  该探测的量测噪声协方差
%     track_id  : 整数  全局唯一轨迹 ID
%     frame_idx : 整数  诞生帧编号
%     cfg       : config_tracker 输出
%
%   OUTPUT
%     track : 轨迹结构体（见下文字段说明）
%
%   初始化策略（Bar-Shalom §3.4.2）：
%     位置初始化为量测值；
%     速度初始化为零（未知），协方差给大不确定性；
%     位置协方差综合量测噪声与先验不确定性。

% 初始状态：速度设为零，待后续帧修正
x0 = [z(1); z(2); cfg.birth_vel_init; cfg.birth_vel_init];

% 初始协方差
pos_var = max(diag(R)) + cfg.P0_pos;   % 量测噪声 + 先验位置不确定性
P0 = diag([pos_var, pos_var, cfg.P0_vel, cfg.P0_vel]);

% 初始化 IMM 模型组
models = imm_init_models(x0, P0, cfg);

%% ── 构建轨迹结构 ─────────────────────────────────────────────────

track.id     = track_id;
track.status = 'tentative';     % 初始状态：暂定
track.age    = 1;               % 年龄（帧数）

% IMM 状态
track.models = models;
track.mu     = cfg.mu0(:);
track.Pi     = cfg.Pi;
track.c_bar_ = cfg.mu0(:);   % 与 imm_mix 输出同维，避免诞生轨迹缺字段

% 融合输出（当前帧即初始值）
track.x = x0;
track.P = P0;

% M/N 历史窗口（长度 = cfg.mn_N，最新在最右）
track.history         = zeros(1, cfg.mn_N);
track.history(end)    = 1;              % 诞生帧算作第一次命中

% 计数器
track.consec_hits   = 1;
track.consec_misses = 0;
track.total_hits    = 1;
track.total_misses  = 0;

% 诞生与最后更新帧
track.birth_frame    = frame_idx;
track.last_hit_frame = frame_idx;

% 轨迹历史（用于最终评估与绘图）
track.state_hist = x0;
track.P_hist     = P0;
track.mu_hist    = cfg.mu0(:);
track.time_hist  = frame_idx;

% PF 可视化日志（保持 struct 数组字段一致）
track.ess_hist      = [];
track.particle_snap = cell(1, 3);
end
