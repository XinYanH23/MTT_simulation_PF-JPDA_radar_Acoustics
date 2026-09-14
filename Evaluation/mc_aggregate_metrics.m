function agg = mc_aggregate_metrics(m_list)
%MC_AGGREGATE_METRICS  将 1×n_mc 的 eval_exp1_metrics 结果聚合成均值/标准差/CI。
if isstruct(m_list) && ~iscell(m_list)
    tmp = cell(1, numel(m_list));
    for i = 1:numel(m_list)
        tmp{i} = m_list(i);
    end
    m_list = tmp;
end
n = numel(m_list);
if n < 1
    error('mc_aggregate_metrics:Empty', 'm_list 为空。');
end

fields = {'ospa_mean', 'ospa_rmse', 'ospa_loc_mean', 'ospa_card_mean', ...
          'ospa_assign_mean', 'card_err_mean', 'peak_confirmed', ...
          'false_confirm_ratio', 'mean_track_life', 'median_track_life', ...
          'fragmentation_count', 'id_switches', 'pos_rmse'};

agg = struct();
agg.n_mc = n;
agg.raw  = struct();
for k = 1:numel(fields)
    f = fields{k};
    x = nan(n, 1);
    for i = 1:n
        if isfield(m_list{i}, f) && ~isempty(m_list{i}.(f))
            x(i) = m_list{i}.(f);
        end
    end
    agg.(f) = mc_scalar_stats(x);
    agg.raw.(f) = x;
end

N = numel(m_list{1}.ospa_k);
ospa = nan(n, N);
ospa_loc = nan(n, N);
ospa_card = nan(n, N);
ospa_asg = nan(n, N);
nest = nan(n, N);
card = nan(n, N);
ntru = m_list{1}.n_true_k;
for i = 1:n
    ospa(i, :) = m_list{i}.ospa_k;
    nest(i, :) = m_list{i}.n_est_k;
    card(i, :) = m_list{i}.card_err_k;
    if isfield(m_list{i}, 'ospa_loc_k')
        ospa_loc(i, :) = m_list{i}.ospa_loc_k;
        ospa_card(i, :) = m_list{i}.ospa_card_k;
        ospa_asg(i, :)  = m_list{i}.ospa_assign_k;
    end
end
agg.ospa_k        = mc_series_stats(ospa);
agg.ospa_loc_k    = mc_series_stats(ospa_loc);
agg.ospa_card_k   = mc_series_stats(ospa_card);
agg.ospa_assign_k = mc_series_stats(ospa_asg);
agg.n_est_k       = mc_series_stats(nest);
agg.card_err_k    = mc_series_stats(card);
agg.n_true_k      = ntru;
end
