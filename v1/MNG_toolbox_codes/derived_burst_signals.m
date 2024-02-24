function derived_burst_signals(app)
%DERIVED_CARDIAC_SIGNALS interpolates the hb_res.dt_instantaneous
%   Detailed explanation goes here
tmp = nan(1,length(app.data));
idx=0;
for i = 1: length(app.data)
    tmp(i) =  app.data(i).tic_multipl;
    if strcmp (app.data(i).name,'BURST INTEGRAL')
        idx = i;
    end
end
if idx == 0
    idx =length(app.data)+1;
end

min_ts = app.data(find(tmp ==1,1)).ts(1);

% app.data(idx).name = "RMS MSNA";
% app.data(idx).unit = "µV";
% app.data(idx).tic_multipl = 1;
% app.data(idx).data = app.burst_res.x';
% app.data(idx).ts = app.burst_res.ts(2:3);
% app.data(idx).derived = true;
%  
% idx = idx+1;
tmp = 0.01:0.01:diff(app.settings.interval(1,:));

app.data(idx).name = 'BURST INTEGRAL';
app.data(idx).unit = 'µV*ms';
app.data(idx).tic_multipl = 0.01/min_ts;
app.data(idx).data = interp1(app.burst_res.t_burst,app.burst_res.burst_int,tmp)';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;

idx = idx+1;

app.data(idx).name = 'BURST DURATION';
app.data(idx).unit = 's';
app.data(idx).tic_multipl = 0.01/min_ts;
app.data(idx).data = interp1(app.burst_res.t_burst,app.burst_res.burst_dur,tmp)';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;

idx = idx+1;

app.data(idx).name = 'BURST amplitude';
app.data(idx).unit = 's';
app.data(idx).tic_multipl = 0.01/min_ts;
app.data(idx).data = interp1(app.burst_res.t_burst,app.burst_res.burst_amp,tmp)';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;

idx = idx+1;

app.data(idx).name = 'BURST INTERBURST INTERVAL';
app.data(idx).unit = 's';
app.data(idx).tic_multipl = 0.01/min_ts;
app.data(idx).data = interp1(app.burst_res.t_burst(1:end-1),app.burst_res.dt_burst,tmp)';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;

update_signal_popup(app)
end

