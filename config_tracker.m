function cfg = config_tracker()
% CONFIG_TRACKER  所有可调参数集中在此，不散落在各模块中。
%
%   cfg = config_tracker()
%
%   参考：
%     Bar-Shalom (2001) Estimation with Applications to Tracking and Navigation
%     Blackman & Popoli (1999) Design and Analysis of Modern Tracking Systems

cfg = struct();

%% ── 时间 ──────────────────────────────────────────────────────────────
cfg.dt = 0.1;               % 采样间隔 (s)

%% ── IMM 模型 ────────────────────────────────────────────────────────
cfg.n_models = 2;

% 模型 1：CV（恒速）— 低机动噪声强度
cfg.cv_q = 0.5;             % 过程噪声强度 (m²/s³)

% 模型 2：CA 近似（高机动）— 高噪声强度
cfg.ca_q = 8.0;

% Markov 转移矩阵  Pi(i,j) = P(切换到模型j | 当前为模型i)
cfg.Pi = [0.95  0.05;
          0.10  0.90];

% 初始模型概率
cfg.mu0 = [0.8; 0.2];

%% ── 初始协方差 ──────────────────────────────────────────────────────
cfg.P0_pos = 50^2;          % 初始位置方差 (m²)
cfg.P0_vel = 20^2;          % 初始速度方差 (m²/s²)

%% ── 测量模型 ────────────────────────────────────────────────────────
cfg.H = [1 0 0 0;
         0 1 0 0];

% 各传感器噪声（若 Detection 中未携带则使用默认值）
cfg.R_radar    = diag([1.0^2, 1.0^2]);
cfg.R_acoustic = diag([8.0^2, 8.0^2]);

%% ── 验证门 ─────────────────────────────────────────────────────────
% χ²(2) 分位数：95% → 5.991，99% → 9.210，99.9% → 13.816
cfg.gate_chi2 = 13.816;     % 选择较宽的门以不遗漏真实目标

%% ── JPDA ────────────────────────────────────────────────────────────
cfg.Pd         = 0.90;      % 检测概率
cfg.lambda_c   = 1e-4;      % 杂波空间密度 (m⁻²)
cfg.max_events = 5000;      % 联合事件枚举上限（超出后降级）

%% ── 轨迹管理 ────────────────────────────────────────────────────────
% M/N 确认逻辑：最近 N 帧中至少 M 次命中则确认
% 原 3/5 过松易假确认；4/6 过严又漏跟。折中 3/5 确认，靠 birth 抑制控冗余。
cfg.mn_M = 3;
cfg.mn_N = 5;

% 删除阈值（连续漏检次数）
cfg.del_tentative = 3;      % Tentative：略放宽，避免刚出生即被删导致反复重建
cfg.del_confirmed = 30;     % Confirmed：3.0 s 保活，给遮挡/声学弱区缓冲

% 声学知情生命周期（默认关，主实验仍用固定 3/5 与 30）
% 开启后：场弱则确认航迹可提前删；场强则 2/5 即可确认。硬上限仍是 del_confirmed。
% 诞生仍只允许未关联雷达点（birth_radar_only），场只做门限，不单独起始。
cfg.ac_informed_delete = false;
cfg.ac_aided_confirm   = false;
cfg.del_confirmed_if_weak = 8;   % 场弱时连续漏检达到该值即可删（<30）
cfg.mn_M_if_strong        = 2;   % 场强时确认所需命中数
cfg.ac_field_weak_eta     = [];  % 空：用 birth_acoustic_eta
cfg.ac_field_strong_eta   = [];  % 空：1e-2（Softmax 真值处约 2e-2）

%% ── 新轨迹诞生 ──────────────────────────────────────────────────────
cfg.birth_vel_init = 0.0;   % 初始速度估计 (m/s)，未知时置零

% 仅用雷达探测诞生新轨迹（true，推荐）。声学点探测密集且定位噪声大，
% 其信息已通过声学概率场（PF 权重）参与跟踪，无需再诞生轨迹，
% 避免每目标周围生成大量并发冗余轨迹。false 时恢复"所有探测均可诞生"。
cfg.birth_radar_only = true;

% 邻近抑制：未关联雷达点若落在任一**已确认**航迹 birth_suppress_radius 内，
% 不新建航迹（视为已有航迹的漏关联/门外交点，而非新目标）。
% 仅对 confirmed 抑制，避免 tentative 杂波封锁真实新目标出生。
% Exp1 碎裂压测：35→50，逼系统保活而非断轨后换 ID 重建。
% 0 关闭。
cfg.birth_suppress_radius = 50;   % m

% 声学真实性门限（论文 3.2.4）：未关联雷达点处声学似然须 > eta_a 才允许起始。
% Softmax 场已归一化；真值处 L 中位数约 2e-2，随机点多数 ≤1e-6。
% 保持 1e-3（勿升高：过严会阻碍真实重建）。设为 0 关闭。
cfg.birth_acoustic_eta = 1e-3;

%% ── 异构传感分离原则 ────────────────────────────────────────────────
% true（推荐）：JPDA 关联与 PF 雷达似然仅使用雷达点迹 Z_k^r，声学信息只经
%               声学概率场（似然场）进入 PF 权重更新，避免同一传感器信息同时
%               进入数据关联与似然加权（double counting）。
% false       ：恢复原行为，雷达+声学点迹混合进入 JPDA。
cfg.jpda_radar_only = true;

%% ── 轨迹合并（并发冗余去重）────────────────────────────────────────
% 每帧对活跃 confirmed 轨迹做一次去重：位置接近且速度方向一致者判为同一
% 目标的重复轨迹，保留命中更多/更年长者，其余标记 deleted。
% 速度门用于避免在目标交叉时误并两个真实目标。
% 注意：birth_vel_init=0 时新生航迹速度估计偏差大，原 merge_vel=5 几乎
% 永不触发合并；放宽至 15，并由 track_merge 对年轻航迹进一步放宽速度门。
cfg.track_merge_enable = true;
cfg.merge_radius       = 35;    % 位置合并阈值 (m)；断轨重建后近邻合并压碎裂
cfg.merge_vel          = 15;    % 速度差合并阈值 (m/s)
cfg.merge_young_age    = 50;    % 年龄 < 该值：仅用位置门（速度估计不可靠）
cfg.merge_tentative    = false; % 勿合并 tentative：虚假 confirmed 会吞掉真实新生航迹

%% ── PDA 量测噪声合成（答辩 B2，默认 false = Bar-Shalom 直加权和）────
% false: R̄ = Σ_d β_d R_d
% true:  R̄ = Σ_d β_d R_d / (1−β_0)   （条件于”已检测”的归一化权）
cfg.pda_Rbar_normalized = false;

%% ── 雷达优先融合（消融实验开关）────────────────────────────────────
% false（默认）：原始融合，所有传感器探测均参与 JPDA/PDA 更新
% true          ：雷达优先融合（Radar-Priority Fusion）
%                 若某轨迹验证门内存在 ≥1 个雷达探测，则本帧仅用雷达探测做 PDA 更新；
%                 若门内无雷达探测，则回退使用声学探测（Acoustic Fallback）。
%                 此开关仅影响 PDA 更新步骤，不改变 JPDA/IMM/EKF/轨迹管理逻辑。
cfg.acoustic_assist_only = false;

%% ── 调试与输出 ───────────────────────────────────────────────────────
cfg.allow_birth    = true;   % false 时未关联探测不诞生新轨迹（验收/调试）
cfg.debug          = false;
cfg.plot_every     = 20;    % 每隔多少帧绘制一次中间结果
cfg.save_mat       = true;

%% -- IMM-PF 粒子滤波参数 -------------------------------------------------
% 主开关：true=使用 PF 替换 KF，false=保持原 KF 流程
cfg.use_pf              = true;

% 每条轨迹每个子模型的粒子数
cfg.pf_N                = 500;

% 重采样阈值：ESS < pf_resample_thresh * N 时触发
cfg.pf_resample_thresh  = 0.5;

% 初始化时从初始 P0 采样粒子的扩散倍率（=1 直接用 P0）
cfg.pf_init_cov_scale   = 1.0;

% PF 更新接口：
%   'legacy_likelihood' 论文主实验：w∝w_*(βλc+ΣβN)*[γL+(1-γ)]，γ=0.3
%   'jpda_posterior'    更严格的贝叶斯接口，见 report/theory_bayes_reformulation/
cfg.pf_update_mode      = 'legacy_likelihood';

% 声学融合权重（仅 legacy_likelihood 使用）
cfg.pf_ac_gamma_mode    = 'fixed';
cfg.pf_ac_gamma_fixed   = 0.3;

% 粒子坐标(m) -> 声学概率场坐标(mm) 的缩放因子
cfg.pf_m2ac_scale       = 1000.0;

% 粒子落在声学概率场网格外时的背景似然
cfg.pf_ac_bg_likelihood = 1e-6;

% 模型似然下限（防止 imm_update_mu 退化）
cfg.pf_lambda_floor     = 1e-300;

end
