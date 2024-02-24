function save_annotate_results(app)

[int_idxs,~] = listdlg('PromptString',{'Please select intervals ',...
    'to be plotted/saved.',''},...
    'SelectionMode','multiple','ListString',app.popup_int_spect.Items);

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg','.eps'});
%path = uigetdir;
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';

int_names = {'full interval'};
borders = nan(length(tmp)+1,2);

tmp =find(vertcat(app.burst_ints.type) == 1);
for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    int_names{i+1} = app.burst_ints(tmp(i)).name;
end

for i = 1:length(int_idxs)
    h = figure('Position', get(0, 'Screensize'),'Visible','off');
    set(h, 'NumberTitle', 'off', ...
    'Name', int_names{int_idxs(i)});

    if ~isnan(app.settings.channel_idx.msna)
        [data_msna,ts_msna,name_msna, unit_msna] = current_signal(app,  app.settings.channel_idx.msna);
        data_msna(:,2) = ts_msna(1):ts_msna(1):ts_msna(2);
        if ~isnan(borders(1,1))
            tmp_borders = [ceil(borders(1)/ts_msna(1)) , floor(borders(2)/ts_msna(1))];
            data_msna = data_msna(tmp_borders(1) : tmp_borders(2),:);
        end
        subplot(4,10,1:8)
        plot(downsample(data_msna(:,2),100), downsample(data_msna(:,1),100),'LineWidth',1.2)
        hold on
        ylabel(['MSNA RAW' ' [' unit_msna ']'])
        xlabel('time [s]')
        title(['MSNA RAW intervall ' int_names{int_idxs(i)}])
        
        subplot(4,10,9:10)
        plot_spec(data_msna(:,1), 1/ts_msna(1), [0.5 0.6470 0.9410], [app.edt_min_freq.Value app.edt_max_freq.Value])
        alpha(0.25)
    end
    
   
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
        subplot(4,10,11:18)
        plot(data_bp(:,2) ,data_bp(:,1) ,'color','black','LineWidth',1)
        hold('on')
        ylabel([name_bp ' [' unit_bp ']'])
        xlabel('time [s]')
        title(['BP interval: ' int_names{int_idxs(i)}]);
        c = [0.4940, 0.1840, 0.5560;...
             0.4660, 0.6740, 0.1880];
        plot(data_bp(ft_idx,2), data_bp(ft_idx,1),'-.','LineWidth',2.5)
        plot(data_bp(sys_idx,2), data_bp(sys_idx,1),'-.','LineWidth',2.5)
        plot(data_bp(ft_idx,2), data_bp(ft_idx,1), '<', 'color', c(1,:), 'markerfacecolor', c(1,:))
        plot(data_bp(sys_idx,2), data_bp(sys_idx,1), '^', 'color', c(2,:), 'markerfacecolor', c(2,:))
        subplot(4,10,19:20)
        plot_spec(data_bp(:,1),1/ts_bp(1),[0.4 0.4 0.4], [app.edt_min_freq.Value app.edt_max_freq.Value])
        alpha(0.25)
    
    end
    
    
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
        subplot(4,10,21:28)
        plot (hb(:,1),hb(:,2),'LineWidth',.7,'Color',[0 0.5 0.2])
        hold on
        %%
        plot(downsample(data_ecg(:,2),20),downsample(data_ecg(:,1),20),'color',[0.7 .2 .3],'LineWidth',1.5)
        hold on
        ylabel([name_ecg ' [' unit_ecg ']'])
        xlabel('time [s]')
        title(['ECG intervall: ' int_names{int_idxs(i)}]);
        plot(data_ecg(rw_idx,2),data_ecg(rw_idx,1), '^', 'color', [0 0.4 0.1], 'markerfacecolor', [0 0.4 0.1])
        hold  off
        subplot(4,10,29:30)
        plot_spec(data_ecg(:,1),1/ts_ecg(1),[1 0.7250 0.6980], [app.edt_min_freq.Value app.edt_max_freq.Value])
        alpha(0.25)
        
    end
    

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
            peak_idx = app.resp_res.idx(:,2);
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
        subplot (4,10,31:38)
        plot (resp(:,1),resp(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2])
    %%%
        hold('on')
    
        plot(data_resp(:,2),data_resp(:,1),'color',[0 .7 .3],'LineWidth',3.5)
        ylabel(['Respiration belt' ' [' unit_resp ']'])
        xlabel('time [s]')
        title(['BELT interval: ' int_names{int_idxs(i)}]);
        plot(data_resp(peak_idx,2),data_resp(peak_idx,1), '^', 'color', [0.5 0.1 1], 'markerfacecolor', [0.5 0.1 1])
        hold('off')
        subplot(4,10,39:40)
        plot_spec(data_resp(:,1),1/ts_resp(1),[0.7250 1 0.6980], [app.edt_min_freq.Value app.edt_max_freq.Value])
        alpha(0.25)
    
    end

    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_annotate_' simple_name(int_names{int_idxs(i)}) '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_annotate_' simple_name(int_names{int_idxs(i)}) '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_annotate_' simple_name(int_names{int_idxs(i)}) '.epsc'])
        end
    end
    close(h)
    app.lbl_working.Text = [num2str(round((100*i)/length(int_idxs))) '% done'];
end

end

function plot_spec(x,fs,PatchCol, xl)

% fs = 1/(xwt_angle(3,1)-xwt_angle(2,1));                                % sample frequency (Hz)
% x=xwt_angle(:,2);
y = fft(x.*hamming(length(x)));
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT

idx = find(f >= 300,1); 
f(idx+1:end) = [];
power(idx+1:end) = [];
plot(f,power,'linewidth',2,'Color','black')
title('Power Using FFT')
xlabel('Frequency Hz')
% ylabel('Power')


patch([0 f f(end)+f(end)-f(end-1)],[-10; power; -10],PatchCol)

xlim(xl)


end