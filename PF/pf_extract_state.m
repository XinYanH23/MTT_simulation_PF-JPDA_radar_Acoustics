function pf = pf_extract_state(pf)
%PF_EXTRACT_STATE  从粒子集计算加权均值和协方差，写入 pf.x / pf.P。
%   供 imm_fuse / imm_update_mu 使用（与 KF model.x / model.P 接口兼容）。

w  = pf.weights(:)';                              % 1 x N
xm = pf.particles * w';                           % nx x 1
dx = pf.particles - xm;                           % nx x N
P  = (dx .* w) * dx';                             % nx x nx
P  = (P + P') / 2;
pf.x = xm;
pf.P = P;
end
