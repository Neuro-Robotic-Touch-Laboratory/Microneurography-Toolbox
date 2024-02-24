function wavelet_analysis(app)
%wavelet_analysis Summary of this function goes here
%   Detailed explanation goes here
cla(app.ax_wavelet_sig1)
cla(app.ax_wavelet_sig2)
cla(app.ax_xwt1)
cla(app.ax_xwt2)
cla(app.ax_xwt3)
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

seriesname={name1, name2};

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
else
    int_idx = int_idx-1;
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = app.burst_ints(tmp(int_idx)).borders;
    borders = [ceil(borders(1)/.01), floor(borders(2)/.01)];
end

d1 = d10(borders(1):borders(2),:);
d2 = d20(borders(1):borders(2),:);

% for i=3
%     if i==1
%         d1=d10(201:floor(length(d10)/2),:);
%         d2=d20(201:floor(length(d10)/2),:);
%     end
%     if i==2
%         d1=d10(floor(length(d10)/2):end,:);
%         d2=d20(floor(length(d10)/2):end,:);
%     end
% 
%     if i==3
%         d1=d10(201:end-201,:);
%         d2=d20(201:end-201,:);
%     end
        



    
   
   

    %% Continuous wavelet transform (CWT)
    % The CWT expands the time series into time
    % frequency space.
    
    % figure('color',[1 1 1])
    %subplot(3,2,3)
    tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];
    wt_ax(app.ax_wavelet_sig1,d1,'S0',yp1,'Dj',1/20,'MaxScale',yp2);
    app.lbl_wavelet_sig1_ttl.Text = seriesname{1};
app.lbl_working.Text = '10 % done';
pause(.05)    
    
    set(app.ax_wavelet_sig1,'xlim',tlim);
    hold(app.ax_wavelet_sig1,'on')
     mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'}); % to change colorscheme
    colormap(app.ax_wavelet_sig1,mycolormap);                                                            %
    colorbar(app.ax_wavelet_sig1,'off')
    ylabel(app.ax_wavelet_sig1,'Period [s]')
app.lbl_working.Text = '20 % done';
pause(.05)
    
    wt_ax(app.ax_wavelet_sig2,d2,'S0',yp1,'Dj',1/20,'MaxScale',yp2);
    app.lbl_wavelet_sig2_ttl.Text = seriesname{2};
app.lbl_working.Text = '30 % done';
pause(.05)    
    set(app.ax_wavelet_sig2,'xlim',tlim)
    hold(app.ax_wavelet_sig2, 'on')
    colormap(app.ax_wavelet_sig2,mycolormap);                                                           % to change colorscheme
    colorbar(app.ax_wavelet_sig2,'off')
    ylabel(app.ax_wavelet_sig2, 'Period [s]')
app.lbl_working.Text = '40 % done';
pause(.05)
    
    %% Cross wavelet transform (XWT)
    % The XWT finds regions in time frequency space where
    % the time series show high common power.
    % assignin('base','d1',d1);
    % assignin('base','d2',d2);
%     try
    % figure('color',[1 1 1])
%     subplot(3,2,2)
%     hold on
    app.lbl_wavelet_xwt1_ttl.Text = ['XWT: ' name1 ' + ' num2str(lag1) ' s ' ' - ' name2 ];
%     title( )
%     hold on
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
    xwt_ax(app.ax_xwt1,D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[100 100],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '50 % done';
pause(.05)
    % D01=[d1(1:end-lag01+1,1) d1(lag01:end,2)];
    % % xwt(d1',d2','ms',16);
    % xwt(D01,d2(1:end-lag01+1,:),'S0',yp1,'ms',yp2,'ArrowDensity',[100 100],'ArrowSize',1,'ArrowHeadSize',1);
%     xlabel('time [s]')
    ylabel(app.ax_xwt1,'Period [s]')
    % hold on
    % colormap(cool)
    % colorMap = jet(256);
    % colormap(colorMap); 
    
    % mycolormap = customcolormap(linspace(0,1,11), {'#860454','#c51b7c','#dc75ab','#f0b7da','#ffdeef','#f8f7f7','#e5f4d9','#b9e084','#7fbc42','#4d921e','#276418'});
    mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
    % colorbar('southoutside');
    colormap(app.ax_xwt1,mycolormap);
    colorbar(app.ax_xwt1)
    
    hold(app.ax_xwt1,'on')
    colorbar(app.ax_xwt1,'off')
    
    xlim(app.ax_xwt1,[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
%     end
app.lbl_working.Text = '60 % done';
pause(.05)
    
    
%     try
%     subplot(3,2,4)
%     hold on
    app.lbl_wavelet_xwt2_ttl.Text = ['XWT: ' name1 ' + ' num2str(lag2) ' s ' ' - ' name2 ];
%     title( )
%     hold on
    lag01=lag2*1/.01;
    D01 = [d1(1:end-lag01+1,1) d1(lag01:end,2)];
    D02 = d2(1:end-lag01+1,:);

    xwt_ax(app.ax_xwt2,D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[50 50],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '70 % done';
pause(.05)
%     xlabel('time [s]')
    ylabel(app.ax_xwt2,'Period [s]')
    % hold on
    % % colormap(turbo)
%     hold on
    colormap(app.ax_xwt2,mycolormap);   % to change colorscheme
    colorbar(app.ax_xwt2,'off')
    
    xlim(app.ax_xwt2,[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
%     end
app.lbl_working.Text = '80 % done';
pause(.05)
%     try
%     subplot(3,2,6)
%     hold on
    app.lbl_wavelet_xwt3_ttl.Text = ['XWT: ' name1 ' + ' num2str(lag3) ' s ' ' - ' name2 ];
%     title( )
%     hold on
    lag01=lag3*1/.01;
    D01 = [d1(1:end-lag01+1,1) d1(lag01:end,2)];
    D02 = d2(1:end-lag01+1,:);
    xwt_ax(app.ax_xwt3,D01,D02,'S0',yp1,'ms',yp2,'ArrowDensity',[50 50],'ArrowSize',1,'ArrowHeadSize',1);
app.lbl_working.Text = '90 % done';
pause(.05)
%     xlabel('time [s]')
    ylabel(app.ax_xwt3,'Period [s]')
    % hold on
    % % colormap(turbo)
%     hold on
    colormap(app.ax_xwt3,mycolormap);       % to change colorscheme
    colorbar(app.ax_xwt3,'off')
    
    xlim(app.ax_xwt3,[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))]);
% end
end