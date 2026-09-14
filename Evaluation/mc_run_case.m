function pack = mc_run_case(allDetections, truth, time, cfg, ac_cfg, opts, seeds)
%MC_RUN_CASE  固定其余参数，按独立种子批量跑跟踪器并汇总指标。
%
%   pack.seeds     1×n_mc
%   pack.per_seed  1×n_mc cell（eval_exp1_metrics 输出）
%   pack.agg       mc_aggregate_metrics 输出

seeds = seeds(:)';
n = numel(seeds);
if ~isfield(opts, 'verbose')
    opts.verbose = false;
end
% 蒙特卡洛默认静默逐帧打印
opts.verbose = false;

m_list = cell(1, n);
t_batch = tic;
for i = 1:n
    fprintf('    [MC %d/%d] rng(%d)\n', i, n, seeds(i));
    rng(seeds(i), 'twister');
    t_one = tic;
    out = run_tracker_once(allDetections, truth, time, cfg, ac_cfg, opts);
    m_list{i} = eval_exp1_metrics(out.tracks, truth, time, out.cfg, ...
                                  out.frame_log, true);
    m_list{i}.seed = seeds(i);
    used = toc(t_batch);
    one  = toc(t_one);
    eta  = (used / i) * (n - i);
    fprintf('         本趟 %.0fs | 已用 %s | 本方案剩余约 %s  (%d/%d)\n', ...
        one, fmt_hms(used), fmt_hms(eta), i, n);
    prog = getenv('TRACKER_PROGRESS_FILE');
    if ~isempty(prog)
        fid = fopen(prog, 'w', 'n', 'UTF-8');
        if fid >= 0
            fprintf(fid, 'MC %d/%d  seed=%d  last=%.0fs  elapsed=%s  eta=%s\n', ...
                i, n, seeds(i), one, fmt_hms(used), fmt_hms(eta));
            fclose(fid);
        end
    end
end

pack = struct();
pack.seeds    = seeds;
pack.per_seed = m_list;
pack.agg      = mc_aggregate_metrics(m_list);
pack.agg.seeds = seeds;
end

function s = fmt_hms(sec)
sec = max(0, round(sec));
h = floor(sec / 3600);
m = floor(mod(sec, 3600) / 60);
s_ = mod(sec, 60);
if h > 0
    s = sprintf('%dh%02dm%02ds', h, m, s_);
else
    s = sprintf('%dm%02ds', m, s_);
end
end
