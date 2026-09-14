function plot_results(tracks, truth, time, metrics, cfg)
% PLOT_RESULTS  跟踪结果可视化（4 子图）。
%
%   plot_results(tracks, truth, time, metrics, cfg)

figure('Color','w','Position',[50 50 1400 900]);
colors_conf = lines(max(1, length(tracks)));

N_truth = length(truth);

%% ── 子图 1：2D 轨迹 ──────────────────────────────────────────────
subplot(2,2,1); hold on; grid on; box on;
title('\bf 2D Trajectories','Interpreter','latex','FontSize',13);
xlabel('X (m)','Interpreter','latex');
ylabel('Y (m)','Interpreter','latex');

% 真值
for j = 1 : N_truth
    raw = truth(j).pos3D;
    plot(raw(:,1), raw(:,2), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
end
h_truth = plot(nan, nan, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ground Truth');

% 估计轨迹（仅显示已确认）
for i = 1 : length(tracks)
    if tracks(i).total_hits < cfg.mn_M, continue; end
    sh = tracks(i).state_hist;
    c  = colors_conf(mod(i-1, size(colors_conf,1))+1, :);
    plot(sh(1,:), sh(2,:), '-', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(sh(1,1), sh(2,1), 'o', 'Color', c, 'MarkerSize', 6, 'HandleVisibility','off');
end
h_est = plot(nan, nan, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Track Estimate');
legend([h_truth, h_est], 'Location', 'best', 'Interpreter','none');
axis equal;

%% ── 子图 2：RMSE 随时间 ─────────────────────────────────────────
subplot(2,2,2); hold on; grid on; box on;
title('\bf Per-Frame RMSE (Confirmed Tracks)','Interpreter','latex','FontSize',13);
xlabel('Time (s)','Interpreter','latex');
ylabel('RMSE (m)','Interpreter','latex');

% 计算逐帧 RMSE（与子图 1 分开重算，便于可视化）
N_frame = length(time);
rmse_ts = nan(1, N_frame);
for k = 1 : N_frame
    errs = [];
    for i = 1 : length(tracks)
        if tracks(i).total_hits < cfg.mn_M, continue; end
        th = tracks(i).time_hist;
        [~,loc] = min(abs(th - k));
        if isempty(loc) || abs(th(loc)-k) > 1, continue; end
        est_xy = tracks(i).state_hist(1:2, loc);
        % 找最近真值
        best_d = inf;
        for j = 1:N_truth
            raw = truth(j).pos3D;
            if k <= size(raw,1)
                d = norm(est_xy - raw(k,1:2)');
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
valid_k = ~isnan(rmse_ts);
if any(valid_k)
    plot(time(valid_k), rmse_ts(valid_k), 'b-', 'LineWidth', 1.5);
end
yline(5, 'r--', '5 m', 'LabelHorizontalAlignment','right');

%% ── 子图 3：轨迹状态统计 ────────────────────────────────────────
subplot(2,2,3); hold on; grid on; box on;
title('\bf Track Status Over Time','Interpreter','latex','FontSize',13);
xlabel('Frame','Interpreter','latex');
ylabel('Active Track Count','Interpreter','latex');

n_confirmed  = zeros(1, N_frame);
n_tentative  = zeros(1, N_frame);

for i = 1 : length(tracks)
    th = tracks(i).time_hist;
    for ki = 1 : length(th)
        k = th(ki);
        if k < 1 || k > N_frame, continue; end
        if tracks(i).total_hits >= cfg.mn_M
            n_confirmed(k) = n_confirmed(k) + 1;
        else
            n_tentative(k) = n_tentative(k) + 1;
        end
    end
end
plot(1:N_frame, n_confirmed, 'b-', 'LineWidth', 2, 'DisplayName', 'Confirmed');
plot(1:N_frame, n_tentative, 'r-', 'LineWidth', 1, 'DisplayName', 'Tentative');
legend('Location','northeast','Interpreter','none');

%% ── 子图 4：模型概率（第一条已确认轨迹）──────────────────────
subplot(2,2,4); hold on; grid on; box on;
title('\bf IMM Model Probabilities (Track \#1)','Interpreter','latex','FontSize',13);
xlabel('Frame','Interpreter','latex');
ylabel('\mu_j','Interpreter','latex');

found = false;
for i = 1 : length(tracks)
    if tracks(i).total_hits >= cfg.mn_M && ~isempty(tracks(i).mu_hist)
        mu_h = tracks(i).mu_hist;   % M×K
        th   = tracks(i).time_hist;
        for j = 1 : size(mu_h, 1)
            lbl = sprintf('Model %d', j);
            plot(th, mu_h(j,:), 'LineWidth', 1.5, 'DisplayName', lbl);
        end
        found = true;
        break
    end
end
if ~found
    text(0.5,0.5,'No confirmed track yet','Units','normalized',...
        'HorizontalAlignment','center');
end
legend('Location','northeast','Interpreter','none');
ylim([0, 1]);

%% ── 输出性能文字 ─────────────────────────────────────────────────
if isstruct(metrics) && isfield(metrics,'pos_rmse') && ~isnan(metrics.pos_rmse)
    annotation('textbox',[0.35 0.01 0.3 0.05],'String',...
        sprintf('Mean Pos RMSE = %.2f m | Confirmed = %d | False = %d', ...
        metrics.pos_rmse, metrics.confirmed_tracks, metrics.false_tracks),...
        'EdgeColor','none','FontSize',10,'HorizontalAlignment','center');
end
end
