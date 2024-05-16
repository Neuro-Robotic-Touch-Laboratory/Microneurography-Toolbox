function derived_bp_signals(app)

sys_idx = nan;
dia_idx = nan;
mea_idx = nan;
tmp = nan(1,length(app.data));
for i = 1: length(app.data)
    tmp(i) =  app.data(i).tic_multipl;
    if strcmp (app.data(i).name,'systolic BP')
        sys_idx = i;
    end
    if strcmp (app.data(i).name, 'diastolic BP')
        dia_idx = i;
    end
    if strcmp (app.data(i).name, 'mean BP')
        mea_idx = i;
    end
end
[data,ts,~, unit] = current_signal(app, app.settings.channel_idx.bldp);
data(:,2) = (ts(1):ts(1):ts(2))';

if isnan(sys_idx)
    sys_idx = length(app.data)+1;
end

max_fs_idx = find(tmp ==1,1);
app.data(sys_idx).name = 'systolic BP';
app.data(sys_idx).unit = unit;
app.data(sys_idx).tic_multipl = 0.005/app.data(max_fs_idx).ts(1);
tmp = 0.005:0.005:diff(app.settings.interval(1,:));

sys_data = data(app.bp_res.systolic_idx,:);
app.data(sys_idx).data = (interp1(sys_data(:,2),sys_data(:,1),tmp))';
tmp_idx = [find(~isnan(app.data(sys_idx).data),1) ,  find(~isnan(app.data(sys_idx).data),1,'last')];
app.data(sys_idx).data(1:tmp_idx(1)-1) = app.data(sys_idx).data(tmp_idx(1));
app.data(sys_idx).data(tmp_idx(2)+1:end) = app.data(sys_idx).data(tmp_idx(2));
app.data(sys_idx).ts = [tmp(1),tmp(end)];
app.data(sys_idx).derived = true;

if isnan(dia_idx)
    dia_idx = length(app.data)+1;
end

app.data(dia_idx).name = 'diastolic BP';
app.data(dia_idx).unit = unit;
app.data(dia_idx).tic_multipl = 0.005/app.data(max_fs_idx).ts(1);
tmp = 0.005:0.005:diff(app.settings.interval(1,:));
dia_data = data(app.bp_res.foot_idx,:);
app.data(dia_idx).data = (interp1(dia_data(:,2),dia_data(:,1),tmp))';
tmp_idx = [find(~isnan(app.data(dia_idx).data),1) ,  find(~isnan(app.data(dia_idx).data),1,'last')];
app.data(dia_idx).data(1:tmp_idx(1)-1) = app.data(dia_idx).data(tmp_idx(1));
app.data(dia_idx).data(tmp_idx(2)+1:end) = app.data(dia_idx).data(tmp_idx(2));
app.data(dia_idx).ts = [tmp(1),tmp(end)];
app.data(dia_idx).derived = true;

if isnan(mea_idx)
    mea_idx = length(app.data)+1;
end

app.data(mea_idx).name = 'mean BP';
app.data(mea_idx).unit = unit;
app.data(mea_idx).tic_multipl = 0.005/app.data(max_fs_idx).ts(1);
app.data(mea_idx).data = mean([app.data(sys_idx).data,app.data(dia_idx).data],2);
app.data(mea_idx).ts = app.data(sys_idx).ts;
app.data(mea_idx).derived = true;

end