%DIAG_TRACK_MAINT  诊断当前 PF 结果：声学似然分布、长航迹匹配
root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(fullfile(root,'Evaluation'));

S = load(fullfile(root,'Step3_IMM_JPDA_PF_Result.mat'), ...
    'tracks', 'truth', 'time', 'cfg', 'ac_pack_hist', 'frame_log');
truth = S.truth; tracks = S.tracks; cfg = S.cfg; N = numel(S.time);

fprintf('=== frame_log ===\n');
fprintf('confirmed mean=%.2f max=%d median=%.1f\n', ...
    mean(S.frame_log.n_confirmed), max(S.frame_log.n_confirmed), median(S.frame_log.n_confirmed));
fprintf('tentative mean=%.2f max=%d\n', ...
    mean(S.frame_log.n_tentative), max(S.frame_log.n_tentative));

%% 声学似然：真值位置 vs 随机点
ac_hist = S.ac_pack_hist;
n_det = 0; n_field = 0;
Lt = []; Lc = [];
rng(1);
for k = 1:N
    pack = ac_hist{k};
    if isempty(pack) || ~isstruct(pack), continue; end
    n_field = n_field + 1;
    if ~(isfield(pack,'detected') && pack.detected), continue; end
    if ~isfield(pack,'L_softmax') || isempty(pack.L_softmax), continue; end
    n_det = n_det + 1;
    for j = 1:numel(truth)
        if k <= size(truth(j).pos3D,1) && ~any(isnan(truth(j).pos3D(k,1:2)))
            z = truth(j).pos3D(k,1:2)';
            Lt(end+1) = query_L(z, pack, cfg); %#ok<AGROW>
        end
    end
    for t = 1:3
        z = 1000 * rand(2,1);
        Lc(end+1) = query_L(z, pack, cfg); %#ok<AGROW>
    end
end
eta = 0;
if isfield(cfg,'birth_acoustic_eta'), eta = cfg.birth_acoustic_eta; end
fprintf('\n=== acoustic field ===\n');
fprintf('ac frames with pack=%d detected=%d\n', n_field, n_det);
if ~isempty(Lt)
    fprintf('L_truth: n=%d median=%.3e p10=%.3e p90=%.3e\n', ...
        numel(Lt), median(Lt), prctile(Lt,10), prctile(Lt,90));
end
if ~isempty(Lc)
    fprintf('L_rand:  n=%d median=%.3e p10=%.3e p90=%.3e\n', ...
        numel(Lc), median(Lc), prctile(Lc,10), prctile(Lc,90));
end
if eta > 0 && ~isempty(Lt) && ~isempty(Lc)
    fprintf('eta_a=%.3e, truth>=eta: %.2f, rand>=eta: %.2f\n', ...
        eta, mean(Lt>=eta), mean(Lc>=eta));
end

%% 各 UAV 最佳匹配确认航迹
fprintf('\n=== per-UAV best track ===\n');
for j = 1:numel(truth)
    best_hits = 0; best_id = -1; best_med = NaN; best_life = NaN; best_birth = NaN;
    for i = 1:numel(tracks)
        if tracks(i).total_hits < cfg.mn_M, continue; end
        th = tracks(i).time_hist; sh = tracks(i).state_hist;
        errs = [];
        for ii = 1:numel(th)
            k = th(ii);
            if k <= size(truth(j).pos3D,1) && ~any(isnan(truth(j).pos3D(k,1:2)))
                errs(end+1) = norm(sh(1:2,ii) - truth(j).pos3D(k,1:2)'); %#ok<AGROW>
            end
        end
        if ~isempty(errs) && median(errs) < 50 && tracks(i).total_hits > best_hits
            best_hits = tracks(i).total_hits;
            best_id = tracks(i).id;
            best_med = median(errs);
            best_life = th(end) - tracks(i).birth_frame + 1;
            best_birth = tracks(i).birth_frame;
        end
    end
    valid = find(~any(isnan(truth(j).pos3D(:,1:2)),2));
    if best_id > 0
        fprintf('UAV%d alive[%d-%d] best T#%d birth=%d hits=%d life=%d med_err=%.1f m\n', ...
            j, valid(1), valid(end), best_id, best_birth, best_hits, best_life, best_med);
    else
        fprintf('UAV%d alive[%d-%d]: NO matching confirmed track\n', j, valid(1), valid(end));
    end
end

%% 短命确认航迹占比
lives = [];
for i = 1:numel(tracks)
    if tracks(i).total_hits >= cfg.mn_M
        lives(end+1) = tracks(i).time_hist(end) - tracks(i).birth_frame + 1; %#ok<AGROW>
    end
end
fprintf('\n=== confirmed lifespan ===\n');
fprintf('n=%d  median=%d  p25=%d  p75=%d  life<20: %.0f%%  life>=100: %.0f%%\n', ...
    numel(lives), median(lives), prctile(lives,25), prctile(lives,75), ...
    100*mean(lives<20), 100*mean(lives>=100));

function La = query_L(z, ac_pack, cfg)
kappa = cfg.pf_m2ac_scale;
bg = cfg.pf_ac_bg_likelihood;
L_field = ac_pack.L_softmax;
try
    F = griddedInterpolant({ac_pack.xc(:), ac_pack.zc(:)}, double(L_field), 'linear', 'none');
    La = F(z(1)*kappa, z(2)*kappa);
    if isnan(La), La = bg; end
catch
    La = bg;
end
La = max(La, bg);
end
