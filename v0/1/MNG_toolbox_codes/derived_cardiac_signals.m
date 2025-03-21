function derived_cardiac_signals(app)
%DERIVED_CARDIAC_SIGNALS interpolates the hb_res.dt_instantaneous
%   Detailed explanation goes here
idx = 0;
tmp = nan(1,length(app.data));
for i = 1: length(app.data)
    tmp(i) =  app.data(i).tic_multipl;
    if strcmp (app.data(i).name,'RR Intervall')
        idx = i;
    end
end

if idx == 0
    idx = length(app.data)+1;
end
max_fs_idx = find(tmp ==1,1);
app.data(idx).name = 'RR Intervall';
app.data(idx).unit = 's';
app.data(idx).tic_multipl = 0.01/app.data(max_fs_idx).ts(1);
tmp = 0.01:0.01:diff(app.settings.interval(1,:));
app.data(idx).data = (interp1(app.hb_res.t_events(1:end-1),app.hb_res.dt_instantaneous,tmp))';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;

idx = idx +1;

app.data(idx).name = 'Heartrate, averaged 10 RR-ints';
app.data(idx).unit = 'BPM';
app.data(idx).tic_multipl = 0.01/app.data(max_fs_idx).ts(1);
tmp = 0.01:0.01:diff(app.settings.interval(1,:));
app.data(idx).data = (interp1(app.hb_res.t_events(1:end-1),HRV.HR(app.hb_res.dt_instantaneous,10),tmp))';
tmp_idx = [find(~isnan(app.data(idx).data),1) ,  find(~isnan(app.data(idx).data),1,'last')];
app.data(idx).data(1:tmp_idx(1)-1) = app.data(idx).data(tmp_idx(1));
app.data(idx).data(tmp_idx(2)+1:end) = app.data(idx).data(tmp_idx(2));
app.data(idx).ts = [tmp(1),tmp(end)];
app.data(idx).derived = true;



update_signal_popup(app)
end

