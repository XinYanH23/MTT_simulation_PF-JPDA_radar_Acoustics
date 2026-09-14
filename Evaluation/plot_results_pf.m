function fig_pf = plot_results_pf(tracks, truth, time, metrics, cfg, ac_pack_hist, assoc_hist)
%PLOT_RESULTS_PF  IMM-PF 系统 8 幅完整可视化。
%
%   plot_results_pf(tracks, truth, time, metrics, cfg)
%   plot_results_pf(tracks, truth, time, metrics, cfg, ac_pack_hist, assoc_hist)
%
%   ac_pack_hist : N×1 cell，每帧的 ac_pack（可选，用于图4）
%   assoc_hist   : N×1 cell，每帧的 assoc  （可选，用于图6）

if nargin < 6; ac_pack_hist = {}; end
if nargin < 7; assoc_hist   = {}; end

N     = length(time);
n_tgt = length(truth);
pf_N  = 500;
if isfield(cfg, 'pf_N'); pf_N = cfg.pf_N; end

COLORS = lines(max(n_tgt + 2, 6));
conf_tracks = select_pf_plot_tracks(tracks, cfg);
n_conf = length(conf_tracks);

fig_pf = figure('Name', 'IMM-PF 跟踪结果', 'Position', [30 30 1500 900], 'Color', 'w');

%% ════════════════════════════════════════════════════════════════════════
%% 图 1：真实轨迹 vs 估计轨迹
%% ════════════════════════════════════════════════════════════════════════
ax1 = subplot(2, 4, 1); hold on; grid on; box on;
title('图1  轨迹对比', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

h_legend = [];
for n = 1:n_tgt
    pos = truth(n).pos3D;
    if isempty(pos); continue; end
    h = plot(pos(:,1), pos(:,2), '-', 'Color', COLORS(n,:), 'LineWidth', 2);
    h_legend(end+1) = h; %#ok
    % 标注起点/终点
    plot(pos(1,1),   pos(1,2),   'o', 'Color', COLORS(n,:), 'MarkerSize', 6, 'MarkerFaceColor', COLORS(n,:));
    plot(pos(end,1), pos(end,2), 's', 'Color', COLORS(n,:), 'MarkerSize', 6, 'MarkerFaceColor', COLORS(n,:));
end
for ti = 1:n_conf
    sh = conf_tracks(ti).state_hist;
    if isempty(sh); continue; end
    c  = COLORS(mod(ti-1, size(COLORS,1))+1, :);
    plot(sh(1,:), sh(2,:), '--', 'Color', c, 'LineWidth', 1.2);
end
% 图例：真值实线 + 估计虚线
plot(NaN, NaN, 'k-',  'LineWidth', 2, 'DisplayName', '真值');
plot(NaN, NaN, 'k--', 'LineWidth', 1.2, 'DisplayName', '估计');
legend('真值', '估计', 'Location', 'best', 'FontSize', 8);

%% ════════════════════════════════════════════════════════════════════════
%% 图 2：IMM 模型概率曲线
%% ════════════════════════════════════════════════════════════════════════
ax2 = subplot(2, 4, 2); hold on; grid on; box on;
title('图2  IMM 模型概率', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('帧'); ylabel('μ');

for ti = 1:min(n_conf, 3)
    mu_h = conf_tracks(ti).mu_hist;   % 2×T
    t_h  = conf_tracks(ti).time_hist;
    if isempty(mu_h) || size(mu_h,1) < 2; continue; end
    c = COLORS(ti,:);
    plot(t_h, mu_h(1,:), '-',  'Color', c, 'LineWidth', 1.8, ...
         'DisplayName', sprintf('M1-T%d', conf_tracks(ti).id));
    plot(t_h, mu_h(2,:), '--', 'Color', c, 'LineWidth', 1.4, ...
         'DisplayName', sprintf('M2-T%d', conf_tracks(ti).id));
    % 机动区间高亮（μ₂ > 0.5）
    mano_mask = mu_h(2,:) > 0.5;
    if any(mano_mask)
        t_mano  = t_h(mano_mask);
        y_lower = zeros(1, sum(mano_mask));
        y_upper = mu_h(2, mano_mask);
        fill([t_mano, fliplr(t_mano)], [y_upper, zeros(1,length(y_upper))], ...
             c, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end
end
yline(0.5, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
ylim([0 1]); legend('Location', 'best', 'FontSize', 7);

%% ════════════════════════════════════════════════════════════════════════
%% 图 3：粒子云演化快照
%% ════════════════════════════════════════════════════════════════════════
ax3 = subplot(2, 4, 3); hold on; grid on; box on;
title('图3  粒子云快照', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

snap_frames = round([N*0.25, N*0.5, N*0.75]);
alphas      = [0.25, 0.45, 0.75];
markers     = {'.', '.', '.'};
sz          = [3, 4, 5];

% 如果 conf_tracks 记录了粒子快照
has_snap = n_conf > 0 && isfield(conf_tracks(1), 'particle_snap') && ...
           ~isempty(conf_tracks(1).particle_snap);
if has_snap
    for si = 1:length(snap_frames)
        for ti = 1:min(n_conf, 3)
            if si <= length(conf_tracks(ti).particle_snap)
                pts = conf_tracks(ti).particle_snap{si};
                if isempty(pts); continue; end
                scatter(pts(1,:), pts(2,:), sz(si), COLORS(ti,:), ...
                    'filled', 'MarkerFaceAlpha', alphas(si));
            end
        end
    end
    text(0.05, 0.93, sprintf('帧 %d / %d / %d', snap_frames), ...
        'Units','normalized', 'FontSize', 8, 'BackgroundColor', [1 1 1 0.6]);
else
    text(0.5, 0.5, {'粒子快照未记录', '请在主循环中', '保存 particle\_snap'}, ...
        'Units','normalized', 'HorizontalAlignment','center', 'FontSize', 9);
end
% 叠加真值
for n = 1:n_tgt
    pos = truth(n).pos3D;
    if isempty(pos); continue; end
    plot(pos(:,1), pos(:,2), 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
end

%% ════════════════════════════════════════════════════════════════════════
%% 图 4：声学概率场热力图
%% ════════════════════════════════════════════════════════════════════════
ax4 = subplot(2, 4, 4); hold on; box on;
title('图4  声学概率场', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('X (mm)'); ylabel('Z (mm)');

ac_last = [];
if ~isempty(ac_pack_hist)
    for ki = N:-1:1
        if ki <= length(ac_pack_hist) && isstruct(ac_pack_hist{ki}) && ...
                isfield(ac_pack_hist{ki}, 'post') && ~isempty(ac_pack_hist{ki}.post)
            ac_last = ac_pack_hist{ki};
            break;
        end
    end
end

if ~isempty(ac_last)
    imagesc(ac_last.xc, ac_last.zc, ac_last.post');
    colormap(ax4, hot); cb = colorbar; cb.Label.String = 'P(x,z)';
    set(ax4, 'YDir', 'normal');
    if isfield(ac_last, 'est_mean') && ~isempty(ac_last.est_mean)
        plot(ac_last.est_mean(1), ac_last.est_mean(2), 'c+', 'MarkerSize', 10, 'LineWidth', 2);
    end
else
    text(0.5, 0.5, '需传入 ac\_pack\_hist 参数', 'Units','normalized', ...
         'HorizontalAlignment','center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
end

% 标注麦克风节点
try
    ac_cfg_v = acoustic_config();
    mic = ac_cfg_v.mic_coords;
    plot(mic(:,1), mic(:,3), 'b^', 'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
         'DisplayName', '麦克风节点');
    legend('Location', 'best', 'FontSize', 7);
catch; end

%% ════════════════════════════════════════════════════════════════════════
%% 图 5：ESS 曲线
%% ════════════════════════════════════════════════════════════════════════
ax5 = subplot(2, 4, 5); hold on; grid on; box on;
title('图5  ESS 曲线', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('帧'); ylabel('ESS');

has_ess = false;
for ti = 1:min(n_conf, 4)
    if isfield(conf_tracks(ti), 'ess_hist') && ~isempty(conf_tracks(ti).ess_hist)
        t_h = conf_tracks(ti).time_hist;
        ess = conf_tracks(ti).ess_hist;
        % 对齐长度
        min_len = min(length(t_h), length(ess));
        plot(t_h(1:min_len), ess(1:min_len), '-', 'Color', COLORS(ti,:), ...
             'LineWidth', 1.4, 'DisplayName', sprintf('T%d', conf_tracks(ti).id));
        has_ess = true;
    end
end
if ~has_ess
    text(0.5, 0.5, 'ESS 未记录（需在主循环保存 ess\_hist）', ...
        'Units','normalized', 'HorizontalAlignment','center', 'FontSize', 8);
end
h_thresh = yline(0.5 * pf_N, 'r--', 'LineWidth', 1.5);
h_thresh.DisplayName = sprintf('阈值 %.0f', 0.5*pf_N);
ylim([0, pf_N * 1.08]);
legend('Location', 'best', 'FontSize', 7);

%% ════════════════════════════════════════════════════════════════════════
%% 图 6：JPDA 关联概率矩阵
%% ════════════════════════════════════════════════════════════════════════
ax6 = subplot(2, 4, 6); box on;
title('图6  JPDA β矩阵', 'FontSize', 10, 'FontWeight', 'bold');

if ~isempty(assoc_hist)
    best_k = 1; best_nd = 0;
    for ki = 1:min(length(assoc_hist), N)
        if ~isempty(assoc_hist{ki}) && isfield(assoc_hist{ki}, 'beta') && ...
                ~isempty(assoc_hist{ki}.beta)
            nd = size(assoc_hist{ki}.beta, 2);
            if nd > best_nd; best_nd = nd; best_k = ki; end
        end
    end
    a = assoc_hist{best_k};
    if ~isempty(a) && isfield(a, 'beta') && ~isempty(a.beta)
        beta_disp = [a.beta0(:), a.beta];
        imagesc(ax6, beta_disp);
        colormap(ax6, parula); colorbar(ax6);
        xlabel(ax6, '量测列（0=β₀漏检）'); ylabel(ax6, '轨迹行');
        title(ax6, sprintf('图6  JPDA β矩阵 (帧%d)', best_k), ...
              'FontSize', 10, 'FontWeight', 'bold');
        xticks(ax6, 0:size(beta_disp,2)+1);
    end
else
    hold(ax6, 'on');
    text(0.5, 0.5, '需传入 assoc\_hist 参数', 'Units','normalized', ...
         'HorizontalAlignment','center', 'FontSize', 9, 'Parent', ax6);
end

%% ════════════════════════════════════════════════════════════════════════
%% 图 7：RMSE 曲线
%% ════════════════════════════════════════════════════════════════════════
ax7 = subplot(2, 4, 7); hold on; grid on; box on;
title('图7  RMSE 曲线', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('帧'); ylabel('RMSE (m)');

rp = [];
if isfield(metrics, 'rmse_pos_k') && ~isempty(metrics.rmse_pos_k)
    rp = metrics.rmse_pos_k(:)';
elseif isfield(metrics, 'pos_rmse')
    rp = compute_rmse_timeseries_pf(conf_tracks, truth, N);
end
if ~isempty(rp)
    valid_rp = ~isnan(rp);
    if any(valid_rp)
        plot(find(valid_rp), rp(valid_rp), 'b-', 'LineWidth', 1.8, 'DisplayName', 'RMSE_{pos}');
    end
end
if isfield(metrics, 'rmse_vel_k') && ~isempty(metrics.rmse_vel_k)
    rv = metrics.rmse_vel_k;
    valid_rv = ~isnan(rv);
    if any(valid_rv)
        plot(find(valid_rv), rv(valid_rv), 'r--', 'LineWidth', 1.4, 'DisplayName', 'RMSE_{vel}');
    end
end
mean_rmse = NaN;
if isfield(metrics, 'pos_rmse') && ~isempty(metrics.pos_rmse)
    mean_rmse = metrics.pos_rmse;
elseif isfield(metrics, 'rmse_pos') && ~isempty(metrics.rmse_pos)
    mean_rmse = metrics.rmse_pos;
end
if ~isnan(mean_rmse)
    yline(mean_rmse, 'b:', 'LineWidth', 1, 'DisplayName', sprintf('均值 %.2f m', mean_rmse));
end
legend('Location', 'best', 'FontSize', 8);

%% ════════════════════════════════════════════════════════════════════════
%% 图 8：OSPA 曲线
%% ════════════════════════════════════════════════════════════════════════
ax8 = subplot(2, 4, 8); hold on; grid on; box on;
title('图8  OSPA (p=2, c=100m)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('帧'); ylabel('OSPA (m)');

ospa_k = [];
if isfield(metrics, 'ospa') && ~isempty(metrics.ospa)
    ospa_k = metrics.ospa(:)';
end
if isempty(ospa_k)
    ospa_k = compute_ospa_online(conf_tracks, truth, N);
end
if ~isempty(ospa_k)
    plot(1:length(ospa_k), ospa_k, 'm-', 'LineWidth', 1.8, 'DisplayName', 'OSPA');
    % 移动均值
    win = min(20, length(ospa_k));
    ospa_smooth = movmean(ospa_k, win);
    plot(1:length(ospa_smooth), ospa_smooth, 'k-', 'LineWidth', 1.0, ...
         'DisplayName', sprintf('均值(w=%d)', win));
    yline(mean(ospa_k), 'k:', 'LineWidth', 1, 'DisplayName', sprintf('均值 %.1f', mean(ospa_k)));
else
    text(0.5, 0.5, 'OSPA 数据不足（确认轨迹<1 或目标<1）', ...
        'Units','normalized', 'HorizontalAlignment','center', 'FontSize', 8);
end
legend('Location', 'best', 'FontSize', 8);
ylim([0, 110]);

%% ── 全局标题 ──────────────────────────────────────────────────────────────
sgtitle(sprintf('IMM-PF 跟踪结果  |  N=%d 粒子  |  %d 目标  |  %d 帧', ...
    pf_N, n_tgt, N), 'FontSize', 12, 'FontWeight', 'bold');

drawnow;
end

%% ─────────────────── 内部辅助：筛选可绘图轨迹 ────────────────────────────
function plot_tracks = select_pf_plot_tracks(tracks, cfg)
plot_tracks = tracks([]);
hits = [tracks.total_hits];
for i = 1:length(tracks)
    if hits(i) >= cfg.mn_M
        plot_tracks(end+1) = tracks(i); %#ok<AGROW>
    end
end
if ~isempty(plot_tracks)
    [~, ord] = sort([plot_tracks.total_hits], 'descend');
    plot_tracks = plot_tracks(ord);
end
end

%% ─────────────────── 内部辅助：逐帧 RMSE（回退） ─────────────────────────
function rmse_ts = compute_rmse_timeseries_pf(plot_tracks, truth, N)
rmse_ts = nan(1, N);
for k = 1:N
    errs = [];
    for ti = 1:length(plot_tracks)
        th = plot_tracks(ti).time_hist;
        [~, loc] = min(abs(th - k));
        if isempty(loc) || abs(th(loc) - k) > 1, continue; end
        est_xy = plot_tracks(ti).state_hist(1:2, loc);
        best_d = inf;
        for j = 1:length(truth)
            raw = truth(j).pos3D;
            if k <= size(raw, 1) && ~any(isnan(raw(k, 1:2)))
                d = norm(est_xy - raw(k, 1:2)');
                if d < best_d, best_d = d; end
            end
        end
        if best_d < 100
            errs(end+1) = best_d^2; %#ok<AGROW>
        end
    end
    if ~isempty(errs)
        rmse_ts(k) = sqrt(mean(errs));
    end
end
end

%% ─────────────────── 内部辅助：在线 OSPA 计算 ────────────────────────────
function ospa_k = compute_ospa_online(conf_tracks, truth, N)
c_ospa = 100; p_ospa = 2;
ospa_k = zeros(1, N);

for k = 1:N
    est = [];
    for ti = 1:length(conf_tracks)
        th = conf_tracks(ti).time_hist;
        sh = conf_tracks(ti).state_hist;
        if isempty(th) || isempty(sh); continue; end
        [~, idx] = min(abs(th - k));
        if abs(th(idx) - k) <= 1 && ~isempty(idx)
            est = [est, sh(1:2, min(idx, size(sh,2)))]; %#ok
        end
    end
    tru = [];
    for n = 1:length(truth)
        if k <= size(truth(n).pos3D, 1) && ~any(isnan(truth(n).pos3D(k,:)))
            tru = [tru, truth(n).pos3D(k, 1:2)']; %#ok
        end
    end
    if isempty(tru);  ospa_k(k) = c_ospa; continue; end
    if isempty(est);  ospa_k(k) = c_ospa; continue; end

    m = size(est, 2);  n_t = size(tru, 2);  mn = max(m, n_t);

    D = zeros(m, n_t);
    for i = 1:m
        for j = 1:n_t
            D(i,j) = min(norm(est(:,i) - tru(:,j)), c_ospa);
        end
    end
    % 贪婪分配（不引入工具箱依赖）
    asgn = 0;
    D2 = D;
    for s = 1:min(m, n_t)
        [~, bi] = min(D2(:));
        [r, c_] = ind2sub(size(D2), bi);
        asgn = asgn + D2(r, c_)^p_ospa;
        D2(r,:) = Inf; D2(:, c_) = Inf;
    end
    ospa_k(k) = ((asgn + c_ospa^p_ospa * abs(m - n_t)) / mn)^(1/p_ospa);
end
end
