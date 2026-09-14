function K = kernel_eval(dist2, sigma_mm, ktype)
%KERNEL_EVAL 空间核函数（用于 NMS-GMM 多模态似然场）。
%
%   K = KERNEL_EVAL(dist2, sigma_mm, ktype)
%
%   输入：
%     dist2    : 到峰值中心的平方距离 (mm^2)，任意尺寸数组
%     sigma_mm : 核宽度 (mm)
%     ktype    : 'gaussian' | 'epanechnikov' | 'laplacian'
%
%   三种核（docx 2.4.3）：
%     gaussian      解析平滑：exp(-d^2 / (2 sigma^2))
%     epanechnikov  紧支撑、无长尾：max(1 - d^2/sigma^2, 0)
%     laplacian     重尾、抗峰值偏移：exp(-|d| / sigma)

s2 = sigma_mm.^2;
switch lower(ktype)
    case 'gaussian'
        K = exp(-dist2 ./ (2 * s2));
    case 'epanechnikov'
        K = max(1 - dist2 ./ s2, 0);
    case 'laplacian'
        K = exp(-sqrt(dist2) ./ sigma_mm);
    otherwise
        error('kernel_eval:unknownKernel', '未知核类型: %s', ktype);
end

end
