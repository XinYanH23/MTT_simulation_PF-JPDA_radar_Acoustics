%% diagnostic_pf_weight_stats.m
%  权重 / 似然 / ESS 精细诊断脚本
%  ─────────────────────────────────────────────────────────────────────
%  不修改任何现有 .m 文件。
%  调用 PF/pf_jpda_update_debug.m（新建独立函数文件）复现
%  pf_jpda_update.m 全部逻辑，同时捕获以下中间量：
%
%    第1项  更新前权重     max/min(w_prior)
%    第2项  更新后未归一化 max/min(w_prior × ell)
%    第3项  归一化后       max/min(w_norm)
%    第4项  每帧动态范围   max(w_norm)/min(w_norm)
%    第5项  log似然        max/min log(ell_fused)
%    第6项  声学似然       max/min/mean(ell_ac)
%    第7项  雷达似然       max/min/mean(ell_radar)
%    第8项  更新前/后 ESS

clear; clc; close all;
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir,'IMM'));  addpath(fullfile(this_dir,'JPDA'));
addpath(fullfile(this_dir,'EKF')); addpath(fullfile(this_dir,'PF'));
addpath(fullfile(this_dir,'TrackManagement')); addpath(fullfile(this_dir,'Utilities'));
addpath(fullfile(this_dir,'Evaluation')); addpath(fullfile(this_dir,'Scenario'));
addpath(fullfile(this_dir,'matlab'));

%% ── 参数 ──────────────────────────────────────────────────────────────────
cfg        = config_tracker();
cfg.use_pf = true;
ac_cfg     = acoustic_config();
prev_post_ac = [];

s2 = fullfile(this_dir, 'Scenario', 'Step2_HeteroDetections.mat');
s1 = fullfile(this_dir, 'Scenario', 'Step1_Data.mat');
if ~exist(s2,'file') || ~exist(s1,'file')
    error('未找到数据文件，请先运行 run_scenario.m / step2_detection_simulation.m');
end
load(s2, 'allDetections', 'time');
data1 = load(s1);
truth = data1.truth;
N_total = length(time);
cfg.dt  = time(2) - time(1);

N_DIAG = min(50, N_total);
fprintf('===== pf_jpda_update 权重/似然/ESS 诊断（前 %d 帧）=====\n\n', N_DIAG);

%% ── 逐帧统计容器 ─────────────────────────────────────────────────────────
s_wp_max  = nan(N_DIAG,1);   s_wp_min  = nan(N_DIAG,1);  % 1
s_wu_max  = nan(N_DIAG,1);   s_wu_min  = nan(N_DIAG,1);  % 2
s_wn_max  = nan(N_DIAG,1);   s_wn_min  = nan(N_DIAG,1);  % 3
s_le_max  = nan(N_DIAG,1);   s_le_min  = nan(N_DIAG,1);  % 5
s_ac_max  = nan(N_DIAG,1);   s_ac_min  = nan(N_DIAG,1);  s_ac_mean = nan(N_DIAG,1);  % 6
s_rd_max  = nan(N_DIAG,1);   s_rd_min  = nan(N_DIAG,1);  s_rd_mean = nan(N_DIAG,1);  % 7
s_eb      = nan(N_DIAG,1);   s_ea      = nan(N_DIAG,1);  % 8
s_rs      = zeros(N_DIAG,1); % 重采样次数

%% ── 初始化 ────────────────────────────────────────────────────────────────
tracks   = struct([]);
next_id  = 1;
snap_frames = round([N_total*0.25, N_total*0.5, N_total*0.75]);

%% =========================================================================
%%  主循环（前 N_DIAG 帧）
%% =========================================================================
for k = 1:N_DIAG

    %% Step 0：声学概率场
    ac_pack_k = [];
    src_pos_k = [];
    for tn = 1:length(truth)
        if k <= size(truth(tn).pos3D,1) && ~any(isnan(truth(tn).pos3D(k,:)))
            src_pos_k = [src_pos_k; truth(tn).pos3D(k,1:3)]; %#ok
        end
    end
    if ~isempty(src_pos_k)
        try
            [ac_pack_k, prev_post_ac] = acoustic_field_to_tracker( ...
                src_pos_k, ac_cfg, prev_post_ac, 'softmax');
        catch ME
            warning('acoustic frame%d: %s', k, ME.message);
        end
    end

    %% Step 1：IMM 预测（不变）
    for i = 1:length(tracks)
        if strcmp(tracks(i).status,'deleted'), continue; end
        [tracks(i).models, c_bar_i] = imm_pf_mix(tracks(i).models, ...
                                                   tracks(i).mu, tracks(i).Pi);
        tracks(i).c_bar_ = c_bar_i;
        tracks(i).models = imm_pf_predict(tracks(i).models);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).c_bar_);
    end

    %% Step 2&3：量测 + JPDA（不变）
    [Z_k, R_list_k, ~] = detections_to_ZR(allDetections{k}, cfg);
    if isempty(tracks)
        active_mask = false(1,0);
    else
        active_mask = ~strcmp({tracks.status},'deleted');
    end
    active_idx = find(active_mask);
    n_active   = length(active_idx);
    if n_active > 0 && size(Z_k,2) > 0
        assoc = jpda_run(tracks(active_idx), Z_k, R_list_k, cfg);
    else
        assoc.valid_mat  = false(n_active, size(Z_k,2));
        assoc.innov_data = cell(n_active, size(Z_k,2));
        assoc.beta       = zeros(n_active, size(Z_k,2));
        assoc.beta0      = ones(n_active,1);
        assoc.clusters   = {};
    end

    %% Step 4：PDA 更新（用 debug 版替换，不修改原文件）
    acc_wp_max=[]; acc_wp_min=[]; acc_wu_max=[]; acc_wu_min=[];
    acc_wn_max=[]; acc_wn_min=[]; acc_le_max=[]; acc_le_min=[];
    acc_rd_max=[]; acc_rd_min=[]; acc_rd_mn=[];
    acc_ac_max=[]; acc_ac_min=[]; acc_ac_mn=[];
    acc_eb=[]; acc_ea=[]; n_rs=0;

    for ti = 1:n_active
        i = active_idx(ti);

        % 提取本轨迹的关联数据（复现 track_imm_pda_update.m 中的提取逻辑）
        beta_row  = assoc.beta(ti, :);
        beta0_val = assoc.beta0(ti);
        d_valid   = find(assoc.valid_mat(ti, :));
        Dv_i      = length(d_valid);
        Z_valid   = zeros(2, Dv_i);
        R_valid   = zeros(2, 2, Dv_i);
        beta_valid= zeros(1, Dv_i);
        for kd = 1:Dv_i
            d = d_valid(kd);
            Z_valid(:,kd)   = Z_k(:,d);
            R_valid(:,:,kd) = R_list_k(:,:,d);
            beta_valid(kd)  = beta_row(d);
        end

        Lambda_j = zeros(cfg.n_models, 1);

        for j = 1:cfg.n_models
            % 调用 debug 版本（新建文件，不修改原函数）
            [tracks(i).models(j), Lambda_j(j), ~, dbg] = pf_jpda_update_debug( ...
                tracks(i).models(j), cfg.H, ...
                Z_valid, R_valid, beta_valid, beta0_val, cfg.lambda_c, ...
                ac_pack_k, cfg);

            acc_wp_max(end+1) = dbg.w_prior_max;    %#ok
            acc_wp_min(end+1) = dbg.w_prior_min;    %#ok
            acc_wu_max(end+1) = dbg.w_unnorm_max;   %#ok
            acc_wu_min(end+1) = dbg.w_unnorm_min;   %#ok
            acc_wn_max(end+1) = dbg.w_norm_max;     %#ok
            acc_wn_min(end+1) = dbg.w_norm_min;     %#ok
            acc_le_max(end+1) = dbg.log_ell_max;    %#ok
            acc_le_min(end+1) = dbg.log_ell_min;    %#ok
            acc_rd_max(end+1) = dbg.ell_radar_max;  %#ok
            acc_rd_min(end+1) = dbg.ell_radar_min;  %#ok
            acc_rd_mn(end+1)  = dbg.ell_radar_mean; %#ok
            acc_ac_max(end+1) = dbg.ell_ac_max;     %#ok
            acc_ac_min(end+1) = dbg.ell_ac_min;     %#ok
            acc_ac_mn(end+1)  = dbg.ell_ac_mean;    %#ok
            acc_eb(end+1)     = dbg.ess_before;     %#ok
            acc_ea(end+1)     = dbg.ess_after;      %#ok
            n_rs = n_rs + dbg.resampled;
        end

        % IMM 模型概率更新 + 融合（与原主循环一致）
        if isfield(tracks(i),'c_bar_')
            c_bar = tracks(i).c_bar_;
        else
            c_bar = tracks(i).Pi' * tracks(i).mu;
        end
        c_bar = max(c_bar(:), 1e-300);
        tracks(i).mu = imm_update_mu(tracks(i).mu, Lambda_j, c_bar);
        [tracks(i).x, tracks(i).P] = imm_fuse(tracks(i).models, tracks(i).mu);

        % ESS 历史记录
        ess_this = min(acc_ea(max(1,end-cfg.n_models+1):end));
        if ~isfield(tracks(i),'ess_hist') || isempty(tracks(i).ess_hist)
            tracks(i).ess_hist = ess_this;
        else
            tracks(i).ess_hist(end+1) = ess_this; %#ok
        end

        % 粒子快照（与主循环一致）
        if ismember(k, snap_frames)
            si = find(snap_frames==k,1);
            [~,mj] = max(tracks(i).mu);
            tracks(i).particle_snap{si} = tracks(i).models(mj).particles(1:2,:);
        end
    end

    % 聚合本帧统计量
    if ~isempty(acc_wp_max)
        s_wp_max(k)  = max(acc_wp_max);   s_wp_min(k)  = min(acc_wp_min);
        s_wu_max(k)  = max(acc_wu_max);   s_wu_min(k)  = min(acc_wu_min);
        s_wn_max(k)  = max(acc_wn_max);   s_wn_min(k)  = min(acc_wn_min);
        s_le_max(k)  = max(acc_le_max);   s_le_min(k)  = min(acc_le_min);
        s_rd_max(k)  = max(acc_rd_max);   s_rd_min(k)  = min(acc_rd_min);  s_rd_mean(k) = mean(acc_rd_mn);
        s_ac_max(k)  = max(acc_ac_max);   s_ac_min(k)  = min(acc_ac_min);  s_ac_mean(k) = mean(acc_ac_mn);
        s_eb(k)      = mean(acc_eb);      s_ea(k)      = mean(acc_ea);
        s_rs(k)      = n_rs;
    end

    %% Step 5：轨迹管理（不变）
    [tracks, next_id] = track_manage(tracks, assoc, Z_k, R_list_k, ...
                                      k, next_id, cfg);
    if ~isempty(tracks) && isfield(tracks,'c_bar_')
        tracks = rmfield(tracks,'c_bar_');
    end
end  % 主循环

%% =========================================================================
%%  输出报告
%% =========================================================================
SEP  = repmat('─',1,72);
fv   = find(~isnan(s_wp_max))';  % 行向量（fill/plot 需要一致维度）

if isempty(fv)
    fprintf('[警告] 前 %d 帧内无活跃轨迹，无诊断数据。\n', N_DIAG);
    return
end

ratio_vec = s_wn_max(fv) ./ max(s_wn_min(fv), 1e-300);

fprintf('\n%s\n第1项  更新前权重  w_prior\n%s\n', SEP, SEP);
fprintf('%-5s  %-16s  %-16s  %-14s\n','帧','max(w_prior)','min(w_prior)','max/min');
fprintf('%s\n',SEP);
for k = fv(:)'
    mx = s_wp_max(k);  mn = s_wp_min(k);
    fprintf('  %3d  %14.6e  %14.6e  %12.2f\n', k, mx, mn, mx/max(mn,1e-300));
end
fprintf('  全局  max=%.6e  min=%.6e\n', max(s_wp_max(fv)), min(s_wp_min(fv)));

fprintf('\n%s\n第2项  更新后（未归一化）unnorm_w = w_prior × ell_fused\n%s\n', SEP, SEP);
fprintf('%-5s  %-16s  %-16s\n','帧','max(unnorm_w)','min(unnorm_w)');
fprintf('%s\n',SEP);
for k = fv(:)'
    fprintf('  %3d  %14.6e  %14.6e\n', k, s_wu_max(k), s_wu_min(k));
end
fprintf('  全局  max=%.6e  min=%.6e\n', max(s_wu_max(fv)), min(s_wu_min(fv)));

fprintf('\n%s\n第3项  归一化后权重  w_norm\n%s\n', SEP, SEP);
fprintf('%-5s  %-16s  %-16s\n','帧','max(w_norm)','min(w_norm)');
fprintf('%s\n',SEP);
for k = fv(:)'
    fprintf('  %3d  %14.6e  %14.6e\n', k, s_wn_max(k), s_wn_min(k));
end
fprintf('  全局  max=%.6e  min=%.6e\n', max(s_wn_max(fv)), min(s_wn_min(fv)));

fprintf('\n%s\n第4项  每帧动态范围  max(w_norm)/min(w_norm)\n%s\n', SEP, SEP);
fprintf('%-5s  %-18s\n','帧','max/min');
fprintf('%s\n',SEP);
for idx = 1:length(fv)
    fprintf('  %3d  %16.4f\n', fv(idx), ratio_vec(idx));
end
fprintf('  汇总：min=%.2f  mean=%.2f  max=%.2f\n', min(ratio_vec), mean(ratio_vec), max(ratio_vec));

fprintf('\n%s\n第5项  log(ell_fused) 对数似然\n%s\n', SEP, SEP);
fprintf('%-5s  %-16s  %-16s  %-12s\n','帧','max log(ell)','min log(ell)','动态范围(差)');
fprintf('%s\n',SEP);
for k = fv(:)'
    fprintf('  %3d  %14.4f  %14.4f  %10.2f\n', k, s_le_max(k), s_le_min(k), s_le_max(k)-s_le_min(k));
end
fprintf('  全局  max=%.4f  min=%.4f  最大动态范围=%.2f\n', ...
    max(s_le_max(fv)), min(s_le_min(fv)), max(s_le_max(fv)-s_le_min(fv)));

fprintf('\n%s\n第6项  声学似然  ell_ac\n', SEP);
fprintf('       [若 ac_pack.detected=false 则恒为 1.0]\n%s\n', SEP);
fprintf('%-5s  %-14s  %-14s  %-14s\n','帧','max','min','mean');
fprintf('%s\n',SEP);
for k = fv(:)'
    fprintf('  %3d  %12.6f  %12.6f  %12.6f\n', k, s_ac_max(k), s_ac_min(k), s_ac_mean(k));
end
ac_all_max = max(s_ac_max(fv));  ac_all_min = min(s_ac_min(fv));
fprintf('  全局  max=%.6f  min=%.6f  mean=%.6f\n', ac_all_max, ac_all_min, mean(s_ac_mean(fv)));
if abs(ac_all_max-1.0)<1e-9 && abs(ac_all_min-1.0)<1e-9
    fprintf('  ★ 警告：ell_ac 全程 = 1.0\n');
    fprintf('    → 声学场未触发（p_detect_joint < det_threshold=%.2f）\n', ac_cfg.det_threshold);
    fprintf('    → 粒子更新完全依赖雷达似然\n');
end

fprintf('\n%s\n第7项  雷达似然  ell_radar\n%s\n', SEP, SEP);
fprintf('%-5s  %-14s  %-14s  %-14s\n','帧','max','min','mean');
fprintf('%s\n',SEP);
for k = fv(:)'
    fprintf('  %3d  %12.6e  %12.6e  %12.6e\n', k, s_rd_max(k), s_rd_min(k), s_rd_mean(k));
end
fprintf('  全局  max=%.6e  min=%.6e  mean=%.6e\n', ...
    max(s_rd_max(fv)), min(s_rd_min(fv)), mean(s_rd_mean(fv)));

fprintf('\n%s\n第8项  ESS（更新前 vs 归一化后，重采样触发前）\n%s\n', SEP, SEP);
fprintf('%-5s  %-14s  %-14s  %-12s  %-10s\n', ...
    '帧','ESS_before','ESS_after','after/N %%','重采样次数');
fprintf('%s\n',SEP);
for k = fv(:)'
    pct = s_ea(k)/cfg.pf_N*100;
    fprintf('  %3d  %12.2f  %12.2f  %10.1f%%  %6d\n', ...
        k, s_eb(k), s_ea(k), pct, s_rs(k));
end
fprintf('  汇总  ESS_before: mean=%.2f  ESS_after: mean=%.2f  总重采样=%d\n', ...
    mean(s_eb(fv)), mean(s_ea(fv)), sum(s_rs));

%% =========================================================================
%%  绘图（8张子图）
%% =========================================================================
fig = figure('Name','pf_jpda_update 诊断','Position',[50 50 1600 900],'Color','w');

ax1 = subplot(2,4,1); hold on; grid on; box on;
title('①  权重 max（前/未归/归）','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('max(w)');
plot(fv, s_wp_max(fv), 'b-o','MarkerSize',3,'LineWidth',1.2,'DisplayName','w_{prior}');
plot(fv, s_wu_max(fv), 'r-s','MarkerSize',3,'LineWidth',1.2,'DisplayName','unnorm');
plot(fv, s_wn_max(fv), 'g-^','MarkerSize',3,'LineWidth',1.2,'DisplayName','norm');
legend('FontSize',7); set(ax1,'YScale','log'); xlim([1,N_DIAG]);

ax2 = subplot(2,4,2); hold on; grid on; box on;
title('②  权重 min（前/未归/归）','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('min(w)');
plot(fv, max(s_wp_min(fv),1e-300), 'b-o','MarkerSize',3,'LineWidth',1.2,'DisplayName','w_{prior}');
plot(fv, max(s_wu_min(fv),1e-300), 'r-s','MarkerSize',3,'LineWidth',1.2,'DisplayName','unnorm');
plot(fv, max(s_wn_min(fv),1e-300), 'g-^','MarkerSize',3,'LineWidth',1.2,'DisplayName','norm');
legend('FontSize',7); set(ax2,'YScale','log'); xlim([1,N_DIAG]);

ax3 = subplot(2,4,3); hold on; grid on; box on;
title('③④  权重动态范围  max/min (w_{norm})','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('max/min');
plot(fv, ratio_vec, 'k-o','MarkerSize',3,'LineWidth',1.5);
yline(1e3, 'b:','LineWidth',1,'DisplayName','1e3');
yline(1e6, 'r--','LineWidth',1,'DisplayName','1e6');
yline(1e9, 'm--','LineWidth',1,'DisplayName','1e9');
legend('FontSize',7,'Location','best'); set(ax3,'YScale','log'); xlim([1,N_DIAG]);
text(0.02,0.97,sprintf('max=%.2e  mean=%.2e',max(ratio_vec),mean(ratio_vec)), ...
    'Units','normalized','FontSize',7,'VerticalAlignment','top','BackgroundColor',[1 1 1 0.7]);

ax4 = subplot(2,4,4); hold on; grid on; box on;
title('⑤  log(ell_{fused}) 范围','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('log(ell)');
plot(fv, s_le_max(fv),'b-','LineWidth',1.4,'DisplayName','max');
plot(fv, s_le_min(fv),'r-','LineWidth',1.4,'DisplayName','min');
fill([fv fliplr(fv)],[s_le_max(fv) fliplr(s_le_min(fv))], ...
     [0.8 0.8 1],'FaceAlpha',0.3,'EdgeColor','none');
legend('FontSize',7); xlim([1,N_DIAG]);

ax5 = subplot(2,4,5); hold on; grid on; box on;
title('⑦  雷达似然  ell_{radar}','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('ell_{radar}');
plot(fv, s_rd_max(fv), 'b-','LineWidth',1.4,'DisplayName','max');
plot(fv, s_rd_min(fv), 'r-','LineWidth',1.4,'DisplayName','min');
plot(fv, s_rd_mean(fv),'g-','LineWidth',1.2,'DisplayName','mean');
legend('FontSize',7); set(ax5,'YScale','log'); xlim([1,N_DIAG]);

ax6 = subplot(2,4,6); hold on; grid on; box on;
title('⑥  声学似然  ell_{ac}','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('ell_{ac}');
plot(fv, s_ac_max(fv), 'b-','LineWidth',1.4,'DisplayName','max');
plot(fv, s_ac_min(fv), 'r-','LineWidth',1.4,'DisplayName','min');
plot(fv, s_ac_mean(fv),'g-','LineWidth',1.2,'DisplayName','mean');
yline(1.0,'k:','LineWidth',2,'HandleVisibility','off');
legend('FontSize',7); xlim([1,N_DIAG]);

ax7 = subplot(2,4,7); hold on; grid on; box on;
title('⑧  ESS：更新前 vs 更新后','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('ESS');
plot(fv, s_eb(fv),'b-o','MarkerSize',3,'LineWidth',1.4,'DisplayName','ESS_{before}');
plot(fv, s_ea(fv),'r-s','MarkerSize',3,'LineWidth',1.4,'DisplayName','ESS_{after}');
yline(cfg.pf_resample_thresh*cfg.pf_N,'k--','LineWidth',1.5, ...
    'DisplayName',sprintf('阈值 %.0f',cfg.pf_resample_thresh*cfg.pf_N));
yline(cfg.pf_N,'k:','LineWidth',1,'DisplayName',sprintf('N=%d',cfg.pf_N));
ylim([0,cfg.pf_N*1.1]); legend('FontSize',7,'Location','best'); xlim([1,N_DIAG]);

ax8 = subplot(2,4,8); hold on; grid on; box on;
title('重采样次数（/帧）','FontSize',9,'FontWeight','bold');
xlabel('帧'); ylabel('次数');
bar(fv, s_rs(fv),'FaceColor',[0.4 0.6 0.9],'EdgeColor','none');
total_rs = sum(s_rs);
text(0.05,0.92,sprintf('总计 %d 次（前%d帧）',total_rs,N_DIAG), ...
    'Units','normalized','FontSize',9,'BackgroundColor',[1 1 1 0.7]);
xlim([0,N_DIAG+1]);

sgtitle(sprintf('pf_jpda_update 诊断  N=%d 粒子  前%d帧', cfg.pf_N, N_DIAG), ...
    'FontSize',12,'FontWeight','bold');
drawnow;

out_dir = fullfile(this_dir,'report');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
saveas(fig, fullfile(out_dir,'pf_weight_stats.png'));
fprintf('\n图像已保存: %s\n', fullfile(out_dir,'pf_weight_stats.png'));
fprintf('===== 诊断完成 =====\n');
