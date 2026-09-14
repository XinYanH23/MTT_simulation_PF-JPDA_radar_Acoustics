function conf_tracks = select_pf_plot_tracks(tracks, cfg)
%SELECT_PF_PLOT_TRACKS  筛选用于可视化的轨迹（已确认或曾确认的轨迹）。
%
%   conf_tracks = select_pf_plot_tracks(tracks, cfg)
%
%   返回：status=='confirmed' 的轨迹，以及已删除但命中次数达到确认门限的轨迹。
%   当 tracks 为空时安全返回空集。

if isempty(tracks)
    conf_tracks = tracks;
    return
end

mask = strcmp({tracks.status}, 'confirmed');

% 已删除但曾经被确认（total_hits >= mn_M）的轨迹也纳入绘图
if isfield(cfg, 'mn_M') && isfield(tracks, 'total_hits')
    deleted_mask = strcmp({tracks.status}, 'deleted');
    hits = [tracks.total_hits];
    mask = mask | (deleted_mask & hits >= cfg.mn_M);
end

conf_tracks = tracks(mask);
end
