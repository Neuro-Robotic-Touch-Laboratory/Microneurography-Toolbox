function update_annotate_axis(app)

int_idx = find(strcmp(app.popup_int_annotate.Value,app.popup_int_annotate.Items));
int_name = app.popup_int_annotate.Value;
if int_idx == 1
    borders = [nan nan];
else
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = app.burst_ints(tmp(int_idx-1)).borders;
end    
xl =[];

cla(app.ax_annotate_1_1)
cla(app.ax_annotate_1_2)
if ~isnan(app.settings.channel_idx.msna)
    [data_msna,ts_msna,name_msna, unit_msna] = current_signal(app,  app.settings.channel_idx.msna);
    data_msna(:,2) = ts_msna(1):ts_msna(1):ts_msna(2);
    if ~isnan(borders(1,1))
        tmp_borders = [ceil(borders(1)/ts_msna(1)) , floor(borders(2)/ts_msna(1))];
        data_msna = data_msna(tmp_borders(1) : tmp_borders(2),:);
    end
    plot(app.ax_annotate_1_1, downsample(data_msna(:,2),10), downsample(data_msna(:,1),10),'LineWidth',1.2)

    app.lbl_annotate_1_1.Text = ['MNG RAW intervall ' int_name];
    plot_spec_ax(app.ax_annotate_1_2, data_msna(:,1), 1/ts_msna(1), [0.5 0.6470 0.9410], [app.edt_min_freq.Value app.edt_max_freq.Value])
    alpha(app.ax_annotate_1_2,0.25)
    xl(end+1,:) = [data_msna(1,2), data_msna(end,2)];
else
    app.lbl_annotate_1_1.Text = 'MNG not selected';
end

cla(app.ax_annotate_2_1)
cla(app.ax_annotate_2_2)
if ~isnan(app.settings.channel_idx.bldp)
    [data_bp,ts_bp,name_bp, unit_bp] = current_signal(app,  app.settings.channel_idx.bldp);
    data_bp(:,2) = ts_bp(1):ts_bp(1):ts_bp(2);
    if ~isnan(borders(1,1))
        tmp_borders = [ceil(borders(1)/ts_bp(1)), floor(borders(2)/ts_bp(1))];
        data_bp = data_bp(tmp_borders(1) : tmp_borders(2), :);
        
        ft_idx = app.bp_res.foot_idx((app.bp_res.foot_idx >= tmp_borders(1)) & (app.bp_res.foot_idx <= tmp_borders(2)));
        ft_idx = ft_idx -tmp_borders(1) +1;
        sys_idx = app.bp_res.systolic_idx((app.bp_res.systolic_idx >= tmp_borders(1)) & (app.bp_res.systolic_idx <= tmp_borders(2)));
        sys_idx = sys_idx -tmp_borders(1) +1; 
        ntch_idx = app.bp_res.notch_idx((app.bp_res.notch_idx >= tmp_borders(1)) & (app.bp_res.notch_idx <= tmp_borders(2)));
        ntch_idx = ntch_idx -tmp_borders(1) +1; 
        dicr_idx = app.bp_res.dicrotic_idx((app.bp_res.dicrotic_idx >= tmp_borders(1)) & (app.bp_res.dicrotic_idx <= tmp_borders(2)));
        dicr_idx = dicr_idx -tmp_borders(1) +1; 
    else
        ft_idx = app.bp_res.foot_idx;
        sys_idx = app.bp_res.systolic_idx;
        ntch_idx = app.bp_res.notch_idx;
        dicr_idx = app.bp_res.dicrotic_idx;
    end

    plot(app.ax_annotate_2_1, data_bp(:,2) ,data_bp(:,1) ,'color','black','LineWidth',1)
    hold(app.ax_annotate_2_1, 'on')

    app.lbl_annotate_2_1.Text = ['BP intervall ' int_name];
    c = [0.4940, 0.1840, 0.5560;...
         0.4660, 0.6740, 0.1880];
    plot(app.ax_annotate_2_1, data_bp(ft_idx,2), data_bp(ft_idx,1),'-.','LineWidth',2.5)
    plot(app.ax_annotate_2_1, data_bp(sys_idx,2), data_bp(sys_idx,1),'-.','LineWidth',2.5)
    plot(app.ax_annotate_2_1, data_bp(ft_idx,2), data_bp(ft_idx,1), '<', 'color', c(1,:), 'markerfacecolor', c(1,:))
    plot(app.ax_annotate_2_1, data_bp(sys_idx,2), data_bp(sys_idx,1), '^', 'color', c(2,:), 'markerfacecolor', c(2,:))

    plot_spec_ax(app.ax_annotate_2_2, data_bp(:,1),1/ts_bp(1),[0.4 0.4 0.4], [app.edt_min_freq.Value app.edt_max_freq.Value])
    alpha(app.ax_annotate_2_2, 0.25)
    xl(end+1,:) = [data_bp(1,2), data_bp(end,2)];
else
    app.lbl_annotate_2_1.Text = 'Blood pressure not selected';
end

cla(app.ax_annotate_3_1)
cla(app.ax_annotate_3_2)

if ~isnan(app.settings.channel_idx.ecg)
    [data_ecg,ts_ecg,name_ecg, unit_ecg] = current_signal(app,  app.settings.channel_idx.ecg);
    data_ecg(:,2) = ts_ecg(1):ts_ecg(1):ts_ecg(2);
    data_ecg(:,1) = lowpass(data_ecg(:,1),40,1/ts_ecg(1));
    if ~isnan(borders(1,1))
        tmp_borders = [ceil(borders(1)/ts_ecg(1)), floor(borders(2)/ts_ecg(1))];
        data_ecg = data_ecg( tmp_borders(1) : tmp_borders(2), :);
        rw_idx = app.hb_res.idx((app.hb_res.idx >= tmp_borders(1)) & (app.hb_res.idx <= tmp_borders(2)));
        rw_idx = rw_idx -tmp_borders(1) +1;
    else
        rw_idx = app.hb_res.idx;
    end
    
 %%
    
    rw_ts = data_ecg(rw_idx,2);
    hb = nan(length(rw_ts)*3,2);
    tmp_idx = 1:3:length(rw_ts)*3;
    hb(tmp_idx,1)= rw_ts;
    hb(tmp_idx+1,1)= rw_ts;
    yr = [min(data_ecg(:,1)) max(data_ecg(:,1))];
    yl = [yr(1)-diff(yr)/20, yr(2)+diff(yr)/20];
    hb(tmp_idx,2)= yl(1);
    hb(tmp_idx+1,2)= yl(2);
    plot (app.ax_annotate_3_1,hb(:,1),hb(:,2),'LineWidth',.7,'Color',[0 0.5 0.2])
    hold (app.ax_annotate_3_1, 'on')
    %%

    
    plot(app.ax_annotate_3_1, downsample(data_ecg(:,2),20),downsample(data_ecg(:,1),20),'color',[0.7 .2 .3],'LineWidth',1.5)

    app.lbl_annotate_3_1.Text = ['ECG intervall ' int_name];

   
    plot(app.ax_annotate_3_1, data_ecg(rw_idx,2),data_ecg(rw_idx,1), '^', 'color', [0 0.4 0.1], 'markerfacecolor', [0 0.4 0.1])
    plot_spec_ax(app.ax_annotate_3_2, data_ecg(:,1),1/ts_ecg(1),[1 0.7250 0.6980], [app.edt_min_freq.Value app.edt_max_freq.Value])
    alpha(app.ax_annotate_3_2, 0.25)
    hold(app.ax_annotate_3_1,'off')
    xl(end+1,:) = [data_ecg(1,2), data_ecg(end,2)];
else
    app.lbl_annotate_1_1.Text = 'ECG not selected';
end

cla(app.ax_annotate_4_1)
cla(app.ax_annotate_4_2)

if ~isnan(app.settings.channel_idx.resp)
    [data_resp,ts_resp,name_resp, unit_resp] = current_signal(app,  app.settings.channel_idx.resp);
    data_resp(:,2) = ts_resp(1):ts_resp(1):ts_resp(2);
    data_resp(:,1) = lowpass(data_resp(:,1),3,1/ts_resp(1));
    if ~isnan(borders(1,1))
        tmp_borders = [ceil(borders(1)/ts_resp(1)), floor(borders(2)/ts_resp(1))];
        data_resp = data_resp(tmp_borders(1) : tmp_borders(2),:);
        peak_ts = app.resp_res.ts((app.resp_res.idx(:,2) >= tmp_borders(1)) & (app.resp_res.idx(:,2) <= tmp_borders(2)));
        peak_idx = app.resp_res.idx((app.resp_res.idx(:,2) >= tmp_borders(1)) & (app.resp_res.idx(:,2) <= tmp_borders(2)),2);
        peak_idx = peak_idx -tmp_borders(1) +1;
    else
        peak_ts = app.resp_res.ts;
        peak_idx = app.resp_res.idx;
    end
    %%%
    
    resp = nan(length(peak_ts)*3,2);
    tmp_idx = 1:3:length(peak_ts)*3;
    resp(tmp_idx,1)= peak_ts;
    resp(tmp_idx+1,1)= peak_ts;
    yr = [min(data_resp(:,1)) max(data_resp(:,1))];

    yl = [yr(1)-diff(yr)/20, yr(2)+diff(yr)/20];
    resp(tmp_idx,2)= yl(1);
    resp(tmp_idx+1,2)= yl(2);
    plot (app.ax_annotate_4_1, resp(:,1),resp(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2])
%%%
    hold(app.ax_annotate_4_1,'on')

    plot(app.ax_annotate_4_1, data_resp(:,2),data_resp(:,1),'color',[0 .7 .3],'LineWidth',3.5)

    app.lbl_annotate_4_1.Text = ['Respiration Belt intervall ' int_name];
    plot(app.ax_annotate_4_1, data_resp(peak_idx,2),data_resp(peak_idx,1), '^', 'color', [0.5 0.1 1], 'markerfacecolor', [0.5 0.1 1])
    hold(app.ax_annotate_4_1,'off')
    plot_spec_ax(app.ax_annotate_4_2, data_resp(:,1),1/ts_resp(1),[0.7250 1 0.6980], [app.edt_min_freq.Value app.edt_max_freq.Value])
    alpha(app.ax_annotate_4_2, 0.25)
    xl(end+1,:) = [data_resp(1,2), data_resp(end,2)];
else
    app.lbl_annotate_1_1.Text = 'Respiration signal not selected';
end
xl = [max(xl(:,1)), min(xl(:,2))];
xlim(app.ax_annotate_1_1, xl)
app.edt_min_time.Value = xl(1);
app.edt_max_time.Value = xl(2);
end



function plot_spec_ax(ax,x,fs,PatchCol, xl)

% fs = 1/(xwt_angle(3,1)-xwt_angle(2,1));                                % sample frequency (Hz)
% x=xwt_angle(:,2);
y = fft(x.*hamming(length(x)));
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT

idx = find(f >= 300,1); 
f(idx+1:end) = [];
power(idx+1:end) = [];
plot(ax,f,power,'linewidth',2,'Color','black')
title(ax,'Power Using FFT')
xlabel(ax,'Frequency Hz')
% ylabel('Power')


patch(ax,[0 f f(end)+f(end)-f(end-1)],[-10; power; -10],PatchCol)

xlim(ax,xl)


end