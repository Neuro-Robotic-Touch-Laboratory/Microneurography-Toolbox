function wavelet_report(app)
%WAVELET_REPORT Summary of this function goes here
%   Detailed explanation goes here
cla(app.ax_wavelet_sig1)
cla(app.ax_wavelet_sig2)
cla(app.ax_xwt1)
cla(app.ax_xwt2)
cla(app.ax_xwt3)

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg'}); % ,'ListString',{'.fig','.jpg','.eps'}
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';


channel_idx = find(strcmp(app.popup_signal1.Value,app.popup_signal_spect.Items));
[data,ts,name1, unit1] = current_signal(app, channel_idx);
data1(:,1) = ts(1):ts(1):ts(2);
data1(:,2) = data;
name1 = char(string(name1));

channel_idx = find(strcmp(app.popup_signal2.Value,app.popup_signal_spect.Items));
[data,ts,name2, unit2] = current_signal(app, channel_idx);
data2(:,1) = ts(1):ts(1):ts(2);
data2(:,2) = data;
name2 = char(string(name2));

if data1(1,1) ~= data2(1,1)
    if data1(1,1) > data2(1,1)
        tmp = find(data2(:,1) == data1(1,1),1);
        data2(1:tmp-1,:) = [];
    else
        tmp = find(data1(:,1) == data2(1,1),1);
        data1(1:tmp-1,:) = [];
    end
end


yp1 = app.edt_min_period.Value;
yp2 = app.edt_max_period.Value;
lag1 = app.edt_lag1.Value;
lag2 = app.edt_lag2.Value;
lag3 = app.edt_lag3.Value;

if app.chkbx_detrend.Value
    data1(:,2) = detrend(data1(:,2));
    data2(:,2) = detrend(data2(:,2));
end

dsf = [int16(.01 / mean(diff(data1(:,1)))), int16(.01 / mean(diff(data2(:,1)))) ];

d01 = [downsample(data1(:,1),dsf(1)),downsample(data1(:,2),dsf(1))];
d02 = [downsample(data2(:,1),dsf(2)),downsample(data2(:,2),dsf(2))];
tmp_ts = 0 : .01 : max([d01(end,1) d02(end,1)]);
d10 = tmp_ts';
d10(:,2) = interp1(d01(:,1),d01(:,2),tmp_ts');
d20 = tmp_ts';
d20(:,2) = interp1(d02(:,1),d02(:,2),tmp_ts');

d10(1,:) = [];
d20(1,:) = [];

clear data data1 data2 d01 d02
int_idx = find(strcmp(app.popup_int_wavelet.Value,app.popup_int_spect.Items));

if int_idx == 1
    borders = [1 size(d10,1)];
    int_name = {'full interval'};
else
    int_idx = int_idx-1;
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = app.burst_ints(tmp(int_idx)).borders;
    borders = [ceil(borders(1)/.01), floor(borders(2)/.01)];
    int_name = {app.burst_ints(tmp(int_idx)).name};
end

d1 = d10(borders(1):borders(2),:);
d2 = d20(borders(1):borders(2),:);

f1 = figure('Position', get(0, 'Screensize'),'Visible','off');

subplot(3,2,1)
yyaxis('left')
plot(d1(:,1),d1(:,2), 'LineWidth',1.5)
ylabel(unit1)
yyaxis('right')
plot (d2(:,1),d2(:,2))
legend (name1, name2)
title ([name1 ' + ' name2])
xlim ([min([d1(1,1),d2(1,1)]), max([d1(end,1), d2(end,1)])])
xlabel('time [s]')
ylabel(unit2)
app.lbl_working.Text = '5 % done';
pause(.05)
    %% Continuous wavelet transform (CWT)
subplot(3,2,3)
tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];
wt(d1,'S0',yp1,'Dj',1/20,'MaxScale',yp2);
app.lbl_working.Text = '10 % done';
pause(.05)
title (name1)
xlim (tlim)    
hold on
mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'}); % to change colorscheme
colormap(mycolormap);                                                            %
colorbar off
ylabel('Period [s]')
app.lbl_working.Text = '20 % done';
pause(.05)

subplot(3,2,5)    
wt(d2,'S0',yp1,'Dj',1/20,'MaxScale',yp2);
app.lbl_working.Text = '30 % done';
pause(.05)
title (name2)
xlim(tlim)
hold on
colormap(mycolormap);                                                           % to change colorscheme
colorbar off
ylabel('Period [s]')
app.lbl_working.Text = '40 % done';
pause(.05)
    
    %% Cross wavelet transform (XWT)
    % The XWT finds regions in time frequency space where
    % the time series show high common power.

subplot(3,2,2)

lag01=lag1*1/.01;
    
if lag01 == 0
D01 = [d1(1:end,1) d1(1:end,2)];
D02 = d2(1:end,:);
% xwt(d1',d2','ms',16);
else
D01 = [d1(1:end-lag01+1,1) d1(lag01:end,2)];
D02 = d2(1:end-lag01+1,:);
% xwt(d1',d2','ms',16);
end
xwt(D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[100 100],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '50 % done';
pause(.05)
    % D01=[d1(1:end-lag01+1,1) d1(lag01:end,2)];
    % % xwt(d1',d2','ms',16);
    % xwt(D01,d2(1:end-lag01+1,:),'S0',yp1,'ms',yp2,'ArrowDensity',[100 100],'ArrowSize',1,'ArrowHeadSize',1);
xlabel('time [s]')
ylabel('Period [s]')
    % hold on
    % colormap(cool)
    % colorMap = jet(256);
    % colormap(colorMap); 
    
    % mycolormap = customcolormap(linspace(0,1,11), {'#860454','#c51b7c','#dc75ab','#f0b7da','#ffdeef','#f8f7f7','#e5f4d9','#b9e084','#7fbc42','#4d921e','#276418'});
mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
    % colorbar('southoutside');
colormap(mycolormap);
colorbar
    
hold on
colorbar off
    
xlim([min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
title(['XWT: ' name1 ' + ' num2str(lag1) ' s ' ' - ' name2 ] )
%     end
    
app.lbl_working.Text = '60 % done';
pause(.05)    
%     try
subplot(3,2,4)
%     hold on
% app.lbl_wavelet_xwt2_ttl.Text = ;

%     hold on
lag01=lag2*1/.01;
D01 = [d1(1:end-lag01+1,1) d1(lag01:end,2)];
D02 = d2(1:end-lag01+1,:);

xwt(D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[50 50],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '70 % done';
pause(.05)
xlabel('time [s]')
ylabel('Period [s]')
    % hold on
    % % colormap(turbo)
%     hold on
colormap(mycolormap);   % to change colorscheme
colorbar off
    
xlim([min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
title(['XWT: ' name1 ' + ' num2str(lag2) ' s ' ' - ' name2 ])
%     end
app.lbl_working.Text = '80 % done'; 
pause(.05)
%     try
subplot(3,2,6)
%     hold on


%     hold on
lag01=lag3*1/.01;
D01 = [d1(1:end-lag01+1,1) d1(lag01:end,2)];
D02 = d2(1:end-lag01+1,:);
xwt(D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[50 50],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '90 % done';
pause(.05)
%     xlabel('time [s]')
ylabel('Period [s]')
    % hold on
    % % colormap(turbo)
%     hold on
colormap(mycolormap);       % to change colorscheme
colorbar off

xlim([min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
title(['XWT: ' name1 ' + ' num2str(lag3) ' s ' ' - ' name2])
app.lbl_working.Text = 'saving plots';
for j = 1 : length(form_idxs)
    switch form_idxs(j)
        case 1
            %f1.Visible = 'on';
            savefig(f1,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_wavelet_analysis_' simple_name(int_name{1,1}) '_SIG_' name1 '-' name2 '.fig'],'compact')
            %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_wavelet_analysis_' simple_name(int_name{1,1}) '_SIG_' name1 '-' name2 '.fig'])
        case 2
            saveas(f1,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_wavelet_analysis_' simple_name(int_name{1,1}) '_SIG_' name1 '-' name2 '.jpeg'])
%         case 3
%             saveas(f1,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_wavelet_analysis_' simple_name(int_name{1,1}) '_SIG_' name1 '-' name2 '.epsc'])
    end
end
close(f1)
end

