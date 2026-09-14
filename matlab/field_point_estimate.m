function [est_mean, est_map] = field_point_estimate(P, XX, ZZ)
%FIELD_POINT_ESTIMATE 从概率场提取点坐标：后验均值与最大后验（MAP）。
%
%   [est_mean, est_map] = FIELD_POINT_ESTIMATE(P, XX, ZZ)
%
%   est_mean = sum_u u * P(u)        一阶矩（重心），平滑、抗抖动
%   est_map  = argmax_u P(u)         最大后验，锐利但易跳峰
%
%   输出均为 1x2 = [x_mm, z_mm]。

P = P / sum(P(:));

est_mean = [sum(XX(:) .* P(:)), sum(ZZ(:) .* P(:))];

[~, ind] = max(P(:));
est_map = [XX(ind), ZZ(ind)];

end
