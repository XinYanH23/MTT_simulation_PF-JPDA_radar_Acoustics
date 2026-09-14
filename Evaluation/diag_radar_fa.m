% quick radar FA stats from Step2
cd(fileparts(fileparts(mfilename('fullpath'))));
load('Scenario/Step2_HeteroDetections.mat','allDetections','time');
load('Scenario/Step1_Data.mat','truth');
N = numel(time); Nt = numel(truth);
n_radar = 0; n_fa = 0; n_tp = 0; frames_r = 0;
for k = 1:N
    dets = allDetections{k};
    if isempty(dets), continue; end
    if ~iscell(dets), dets = num2cell(dets); end
    gt = zeros(0,2);
    for j = 1:Nt
        if k <= size(truth(j).pos3D,1) && ~any(isnan(truth(j).pos3D(k,1:2)))
            gt(end+1,:) = truth(j).pos3D(k,1:2); %#ok<AGROW>
        end
    end
    has_r = false;
    for d = 1:numel(dets)
        od = dets{d};
        if isempty(od), continue; end
        try
            st = char(od.ObjectAttributes.SensorType);
        catch
            continue
        end
        if ~strcmpi(st,'Radar'), continue; end
        has_r = true;
        n_radar = n_radar + 1;
        m = od.Measurement(:);
        xy = m(1:2)';
        if isempty(gt)
            n_fa = n_fa + 1;
        else
            dd = sqrt(sum((gt - xy).^2, 2));
            if min(dd) < 15
                n_tp = n_tp + 1;
            else
                n_fa = n_fa + 1;
            end
        end
    end
    if has_r, frames_r = frames_r + 1; end
end
fprintf('radar dets=%d  TP~%d  FA~%d (15m gate)\n', n_radar, n_tp, n_fa);
fprintf('FA rate among radar dets=%.1f%%\n', 100*n_fa/max(n_radar,1));
fprintf('FA/frame (all %d frames)=%.3f; (radar frames %d)=%.3f\n', ...
    N, n_fa/N, frames_r, n_fa/max(frames_r,1));
fprintf('True dets/frame approx=%.2f\n', n_tp/N);
fprintf('JPDA lambda_c=1e-4 over 1e6 m^2 map => %.0f clutter/frame if whole map\n', 1e-4*1e6);
fprintf('FalseAlarmRate setting in Step2 = 1e-6 (per resolution cell, Toolbox)\n');
