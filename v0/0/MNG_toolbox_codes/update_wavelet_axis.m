function update_wavelet_axis(app)

cla(app.ax_wavelet_sig1)
cla(app.ax_wavelet_sig2)
cla(app.ax_xwt1)
cla(app.ax_xwt2)
cla(app.ax_xwt3)

channel_idx = find(strcmp(app.popup_signal1.Value,app.popup_signal_spect.Items));
[data,ts,name1, unit1] = current_signal(app, channel_idx);
data1(:,1) = ts(1):ts(1):ts(2);
data1(:,2) = data;

channel_idx = find(strcmp(app.popup_signal2.Value,app.popup_signal_spect.Items));
[data,ts,name2, unit2] = current_signal(app, channel_idx);
data2(:,1) = (ts(1):ts(1):ts(2));
data2(:,2) = data;

if app.chkbx_detrend.Value
    data1(:,1) = detrend(data1(:,1));
    data2(:,1) = detrend(data2(:,1));
end

dsf = [int16(.01 / data1(1,1)), int16(.01 / data2(1,1)) ];

d10 = [downsample(data1(:,1),dsf(1)),downsample(data1(:,2),dsf(1))];
d20 = [downsample(data2(:,1),dsf(2)),downsample(data2(:,2),dsf(2))];
clear data data1 data2
int_idx = find(strcmp(app.popup_int_wavelet.Value,app.popup_int_spect.Items));
if int_idx == 1
    borders = [1 min([size(d10,1),size(d20,1)])];
else
    int_idx = int_idx-1;
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = app.burst_ints(tmp(int_idx)).borders;
    borders = [ceil(borders(1)/.01), floor(borders(2)/.01)];
end

d1 = d10(borders(1):borders(2),:);
d2 = d20(borders(1):borders(2),:);

yyaxis(app.ax_wave_signals,'left')
plot(app.ax_wave_signals,d1(:,1),d1(:,2), 'LineWidth',1.5)

yyaxis(app.ax_wave_signals,'right')
plot (app.ax_wave_signals, d2(:,1),d2(:,2))
% legend (app.ax_wave_signals,{name1, name2})
% title (app.ax_wave_signals, [name1{1,1} ' + ' name2{1,1}])
xlim (app.ax_wave_signals, [min([d1(1,1),d2(1,1)]), max([d1(end,1), d2(end,1)])])
app.lbl_signal1.Text = name1;
app.lbl_signal2.Text = name2;







