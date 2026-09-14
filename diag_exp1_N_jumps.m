%% 把“N 为什么抖”拆成：真值阶梯、单次种子跳变、同一时刻 25 种子的分歧
clear; clc;

this_dir = fileparts(mfilename('fullpath'));
mat_b = fullfile(this_dir, 'methods_v2_bayes', 'experiments', ...
    'exp1_acoustic_ablation', 'exp1_results.mat');
S = load(mat_b, 'mc_B', 'time');
t = S.time(:)';
if max(t) < 1; t = (0:numel(t)-1) * 0.1; end

ps = S.mc_B.per_seed;
n_mc = numel(ps);
N = numel(ps{1}.n_est_k);
nest = zeros(n_mc, N);
ntrue = ps{1}.n_true_k(:)';
for i = 1:n_mc
    nest(i, :) = ps{i}.n_est_k;
end

fprintf('==== 真值 N_true 变化点 ====\n');
dtrue = find(diff([ntrue(1) ntrue]) ~= 0);
for k = dtrue
    fprintf('  t=%.1fs  N_true: %d -> %d\n', t(k), ntrue(max(k-1,1)), ntrue(k));
end

fprintf('\n==== 单次种子：N 跳了多少次（真值只跳 %d 次）====\n', numel(dtrue));
nj = zeros(1, n_mc);
for i = 1:n_mc
    nj(i) = sum(diff(nest(i,:)) ~= 0);
end
fprintf('  每种子跳变次数: min %d  median %d  max %d  (真值 %d)\n', ...
    min(nj), median(nj), max(nj), numel(dtrue));

% 挑一条中位跳变的种子，列出每一次 +1/-1
[~, i0] = min(abs(nj - median(nj)));
x = nest(i0, :);
chg = find(diff(x) ~= 0) + 1;
fprintf('\n==== seed #%d (跳 %d 次) 的每一次整数跳 ====\n', i0, nj(i0));
for k = chg
    fprintf('  t=%5.1fs  N: %d -> %d   (N_true=%d)\n', t(k), x(k-1), x(k), ntrue(k));
end

% 平台期 t=10–30s，真值=5
mask = (t >= 10) & (t <= 30) & (ntrue == 5);
fprintf('\n==== 平台期 (10-30s, N_true=5)：25 种子每帧的 N 分布 ====\n');
plat = nest(:, mask);
fprintf('  报 3 的帧占比: %.1f%%\n', 100*mean(plat(:)==3));
fprintf('  报 4 的帧占比: %.1f%%\n', 100*mean(plat(:)==4));
fprintf('  报 5 的帧占比: %.1f%%\n', 100*mean(plat(:)==5));
fprintf('  报 6+ 的帧占比: %.1f%%\n', 100*mean(plat(:)>=6));
fprintf('  同一帧 25 种子的标准差(均值): %.2f\n', mean(std(nest(:, mask), 0, 1)));

% 一条种子在平台期自己抖多少
fprintf('\n==== 平台期内，单条种子自己是否也在抖 ====\n');
self_j = zeros(1, n_mc);
for i = 1:n_mc
    xi = nest(i, mask);
    self_j(i) = sum(diff(xi) ~= 0);
end
fprintf('  平台期每种子跳变: min %d  median %d  max %d\n', ...
    min(self_j), median(self_j), max(self_j));

fprintf('\n==== 方法 B Field 种子级碎裂 / 虚警 / 寿命 ====\n');
frag = zeros(1, n_mc); fa = zeros(1, n_mc); life = zeros(1, n_mc);
for i = 1:n_mc
    frag(i) = ps{i}.fragmentation_count;
    fa(i)   = ps{i}.false_confirm_ratio;
    life(i) = ps{i}.median_track_life;
end
fprintf('  碎裂次数:  mean %.1f  (min %d max %d)\n', mean(frag), min(frag), max(frag));
fprintf('  虚警确认:  mean %.1f%%\n', 100*mean(fa));
fprintf('  中位寿命:  mean %.1f 帧 = %.1f s\n', mean(life), mean(life)*0.1);
