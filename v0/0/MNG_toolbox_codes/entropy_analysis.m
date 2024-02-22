function entropy_res = entropy_analysis(app)
%entropy_analysis Summary of this function goes here
%   Detailed explanation goes here


channel_idx = find(strcmp(app.popup_entropy_signal.Value,app.popup_signal_spect.Items));
[data,ts,name, unit] = current_signal(app, channel_idx);
data(:,2) = ts(1):ts(1):ts(2);

tmp =find(vertcat(app.burst_ints.type) == 1);
borders(1,:) = [1 size(data,1)];
int_name{1} = 'full';
for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    borders(i+1,:) = [ceil(borders(i+1,1)/ts(1)), floor(borders(i+1,2)/ts(1))];
    int_name{i+1} = app.burst_ints(tmp(i)).name;
end

cfg            = [];
switch app.popup_method.Value
    case app.popup_method.Items{1}
        cfg.method = 'PE'; % compute permutation entropy
    case app.popup_method.Items{2}
        cfg.method = 'opdPE'; % compute permutation entropy
    case app.popup_method.Items{3}
        cfg.method = 'all';  % compute all implemented ordinal-patterns-based measures
    case app.popup_method.Items{4}
        cfg.method = 'CE'; % we compute conditional entropy of ordinal patterns
    case app.popup_method.Items{5}
        cfg.method = 'rePE'; % compute robust permutation entropy
    case app.popup_method.Items{6}
        cfg.method = 'PEeq'; % compute permutation entropy for ordinal patterns with tied ranks
end
if strcmp (cfg.method, 'opdPE')
    cfg.orderSeq = app.edt_order_seq.Value;
end
if strcmp (cfg.method, 'rePE')
    cfg.lowerThreshold = 0.2;    % the distance that is considered negligible between points
    cfg.upperThreshold = 100;    % the distance between points most probably related to artifact
end
if strcmp (cfg.method, 'all')
    cfg.plot       = true;
else
    cfg.plot       = false;
end
cfg.order      = app.edt_order.Value;    % ordinal pattens of order 3 (4-points ordinal patterns)
cfg.delay      = app.edt_delay.Value;    % delay 2 between points in ordinal patterns              % (one point between successive points in ordinal patterns)
cfg.windowSize = app.edt_window_size.Value;  % window size = 512 time steps
% cfg.plot       = false;
cfg.units      = 'seconds';
% cfg.units       = unit;
%cfg.time       = t_msna; % OPTIONAL time axis for plotting
% cfg.units      = 'seconds';         % OPTIONAL units of time for plotting

entropy_res = struct('outdata',[],'cfg',[],'int_name',[],'signal_name',[]);
% entropy_fix_res = struct('outdata',[],'cfg',[],'int_name',[],'signal_name',[]);
for i = 1: size(borders,1)
    cfg.time = data(borders(i,1):borders(i,2),2)';
%     try
    entropy_res(i).outdata = OPanalysis( cfg, data(borders(i,1):borders(i,2),1)' );
%     catch
%         entropy_res(i).outdata = nan;
%     end
    %entropy_res(i).outdata = OPanalysis_fixed( cfg, data(borders(i,1):borders(i,2),1)' );
    entropy_res(i).cfg = cfg;
    entropy_res(i).int_name = int_name{i};
    entropy_res(i).signal_name = name;
    entropy_res(i).channel_idx = channel_idx;
    
%     entropy_fix_res(i).outdata = OPanalysis( cfg, data(borders(i,1):borders(i,2),1)' );
%     %entropy_res(i).outdata = OPanalysis_fixed( cfg, data(borders(i,1):borders(i,2),1)' );
%     entropy_fix_res(i).cfg = cfg;
%     entropy_fix_res(i).int_name = int_name{i};
%     entropy_fix_res(i).signal_name = name;
%     entropy_fix_res(i).channel_idx = channel_idx;
    app.lbl_working.Text = [num2str(round((100*i)/size(borders,1))) '% done'];
    pause(0.01)
end

end



