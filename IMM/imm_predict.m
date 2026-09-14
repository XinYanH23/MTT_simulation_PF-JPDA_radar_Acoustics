function models = imm_predict(models)
% IMM_PREDICT  对每个模型执行 EKF 线性预测步骤。
%
%   models = imm_predict(models)
%
%   INPUT / OUTPUT
%     models : 1×M struct，每元素含 .x .P .F .Q
%              输入为混合后状态；输出为预测后状态。
%
%   预测方程（线性，等价于 KF 预测）：
%
%     x̂_j(−) = F_j · x̄_{0j}
%     P_j(−)  = F_j · P̄_{0j} · F_j' + Q_j

for j = 1 : length(models)
    F = models(j).F;
    Q = models(j).Q;

    models(j).x = F * models(j).x;
    P            = F * models(j).P * F' + Q;
    models(j).P  = (P + P') / 2;     % 强制对称
end
end
