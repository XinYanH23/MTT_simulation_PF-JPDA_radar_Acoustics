function models = imm_init_models(x0, P0, cfg)
% IMM_INIT_MODELS  为一条新轨迹初始化 IMM 模型组。
%
%   models = imm_init_models(x0, P0, cfg)
%
%   INPUTS
%     x0  : 4×1  初始状态 [px py vx vy]'
%     P0  : 4×4  初始协方差
%     cfg : config_tracker 输出结构体
%
%   OUTPUT
%     models : 1×M struct 数组，每元素含 .x .P .F .Q
%
%   运动模型（均为 4 维状态 [px py vx vy]'）
%   ─────────────────────────────────────────
%   Model 1 — CV（Constant Velocity）
%     F = [1  0  dt  0 ;     Q 用离散白噪声加速度模型（低 q）
%          0  1  0   dt;
%          0  0  1   0 ;
%          0  0  0   1 ]
%
%   Model 2 — 高机动 CV（等效 CA 近似）
%     F 相同，Q 用更大的 q，体现加速度扰动
%     物理含义：无人机在机动时加速度可视为
%     零均值高斯白噪声，q 越大，模型越"宽松"
%
%   离散白噪声加速度 Q（Bar-Shalom 2001, eq. 6.3.3-1）：
%
%     Q = q * [dt³/3   0      dt²/2   0   ]
%             [0       dt³/3  0       dt²/2]
%             [dt²/2   0      dt      0   ]
%             [0       dt²/2  0       dt  ]

dt       = cfg.dt;
n_models = cfg.n_models;

% PF 模式时额外字段：.particles .weights .N（use_pf=false 时置空/0）
use_pf = isfield(cfg, 'use_pf') && cfg.use_pf;
pf_N   = 0;
if use_pf;  pf_N = cfg.pf_N;  end

models = struct('x',         cell(1, n_models), ...
                'P',         cell(1, n_models), ...
                'F',         cell(1, n_models), ...
                'Q',         cell(1, n_models), ...
                'particles', cell(1, n_models), ...
                'weights',   cell(1, n_models), ...
                'N',         cell(1, n_models));

% 状态转移矩阵（两个模型相同，区别仅在 Q）
F = [1  0  dt  0 ;
     0  1  0   dt;
     0  0  1   0 ;
     0  0  0   1 ];

q_list = [cfg.cv_q, cfg.ca_q];

for j = 1 : n_models
    q = q_list(j);

    Q = q * [dt^3/3   0        dt^2/2   0     ;
             0        dt^3/3   0        dt^2/2;
             dt^2/2   0        dt       0     ;
             0        dt^2/2   0        dt    ];

    models(j).x = x0;
    models(j).P = P0;
    models(j).F = F;
    models(j).Q = Q;

    if use_pf
        pf_j              = pf_init(x0, P0, pf_N, cfg.pf_init_cov_scale);
        models(j).particles = pf_j.particles;
        models(j).weights   = pf_j.weights;
        models(j).N         = pf_N;
    else
        models(j).particles = [];
        models(j).weights   = [];
        models(j).N         = 0;
    end
end
end
