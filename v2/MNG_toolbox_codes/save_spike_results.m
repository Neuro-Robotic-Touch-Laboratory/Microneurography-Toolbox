function save_spike_results(app)

[int_idxs,~] = listdlg('PromptString',{'Please select intervals ',...
    'to be plotted/saved.',''},...
    'SelectionMode','multiple','ListString',app.popup_int_spike.Items);

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg','.eps'});

answer = questdlg('print rate panel?', ...
	'print rates', ...
	'yes','no','return and set rate panel','yes');

switch answer
    case 'yes'
        print_rate = true;
    case 'no'
        print_rate = false;
    case 'return and set rate panel'
        return
end

answer = questdlg('print lag panel?', ...
	'print lags', ...
	'yes','no','yes');

switch answer
    case 'yes'
        print_lag = true;
    case 'no'
        print_lag = false;
end

%path = uigetdir;
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';

[data,ts,name, unit] = current_signal(app, app.settings.channel_idx.msna);
data(:,2) = (ts(1) : ts(1) : ts(2))';
name = char(string(name));

int_names = app.popup_int_spike.Items;

if print_rate
    printcell = cell([]);
    hp = figure('Visible', 'off');
    printax = gca;
    printax.XTick = 1:length(app.lstbx_ints.Value);
    printax.XTickLabel = app.lstbx_ints.Value;
    printax.XTickLabelRotation = 30;
    title(printax,'firing rates')
    hold(printax,"on")

    for i= 1: length(app.lstbx_ints.Value)
        printcell{1,i+1} = app.lstbx_ints.Value{1,i};
    end

    bin_step = app.edt_step.Value;
    nbins = app.edt_max.Value;
    cols = distinguishable_colors(max(app.spike_res.cluster),[1 1 1; 1 0 0; 0 0 1; 0.80 0.95 0.90]);
    for i= 1:max(app.spike_res.cluster)
        printcell{i+1,1} = ['Cluster ' num2str(i)];

        h = figure('Position', get(0, 'Screensize'),'Visible', 'off');
        
        subplot(1,2,1)
        mrate = nan(length(app.lstbx_ints.Value),1);
        
        for j = 1 :length(app.lstbx_ints.Value)
            tmp_idx = find(strcmp(app.lstbx_ints.Value{j},app.lstbx_ints.Items));
            mrate(j) = app.spike_res.spike(tmp_idx).m_fr(i);
            printcell{i+1,j+1} = app.spike_res.spike(tmp_idx).m_fr(i);
        end
        ax = gca;
        
        plot (ax,mrate, 'LineWidth', 2, 'LineStyle','--', 'Color',cols(i,:),'Marker','diamond','MarkerSize',10,'MarkerFaceColor',cols(i,:))
        plot (printax,mrate, 'LineWidth', 2, 'LineStyle','--', 'Color',cols(i,:),'Marker','diamond','MarkerSize',10,'MarkerFaceColor',cols(i,:))
        ax.XTick = 1:size(mrate,1);
        ax.XTickLabel = app.lstbx_ints.Value;
        ax.XTickLabelRotation = 30;
        title(ax,['Cluster ' num2str(i) ' firing rates'])
        ylabel(ax,'[Hz]')
        subplot (1,2,2)

        times = diff(app.spike_res.spike_ts((app.spike_res.cluster == i) & app.spike_res.use_spikes(:,1)')/1000)*1000;
        multi_isi = nnz(times < 3); 
        [N,X]=hist(times,0:bin_step:nbins);
        ax = gca;
        xlim (ax,'manual')
        bar(ax,X(1:end-1),N(1:end-1))
        xlim(ax,[0 nbins]);
        title(ax,['global ISI ' num2str(multi_isi) ' spikes in < 3ms'])

        for k = 1 : length(form_idxs)
            switch form_idxs(k)
                case 1
                    h.Visible = 'on';
                    savefig(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.fig'],'compact')
                    %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.fig'])
                case 2
                    saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.jpeg'])
                case 3
                    saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.epsc'])
            end
        end
        close(h)
    end
    writecell(printcell,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) 'spike_results.xls'], 'Sheet', 'rate panel')
    legend(printax, printcell(2:end,1))
    tmp = {'D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF','AG','AH','AI','AJ','AK','AL','AM','AN','AO','AP','AQ','AR','AS','AT','AU','AV','AW','AX','AY','AZ'};
    xlswritefig(hp, [path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) 'spike_results.xls'], 'rate panel', [tmp{length(app.lstbx_ints.Value)} '1'])
    close(hp)
    delete('h')
    delete('hp')
    clear printcell
end        

if print_lag
    
    h_all = figure('Position', get(0, 'Screensize'),'Visible', 'off');
    ax_all = gca;
    cla(ax_all)
    hold (ax_all, 'on')

    h_all_n = figure('Position', get(0, 'Screensize'),'Visible', 'off');
    ax_all_n = gca;
    cla(ax_all_n)
    hold (ax_all_n, 'on')

    leg = {};
    edges = -0.5 :0.005:3;
    
    spike_ts = app.spike_res.spike_ts/1000;
    spike_clust = app.spike_res.cluster;
    use_spike = app.spike_res.use_spikes;
    all_hb_ts = app.hb_res.t_events;
    beats_use = app.hb_res.use_beats;
    cols = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    plts = [];
    plts_n = [];
    for i  = 1: length(int_idxs)
        h = figure('Position', get(0, 'Screensize'),'Visible', 'off');
        ax_int = gca;
        int_idx = int_idxs(i);
        
       
        spks = {};
        spike_use = use_spike(:,int_idx);
        hb_ts = all_hb_ts(beats_use(:,int_idx));
        for j = 1 : length(hb_ts) 
            spks{1,j} = spike_ts(spike_ts>=(hb_ts(j)-0.5) & spike_ts<=(hb_ts(j)+3) & spike_use');
            spks{1,j} = spks{1,j}-hb_ts(j);
            
        end
               
        num_beats = app.edt_num_beat.Value;
        mov_mean = app.edt_lag_movmean.Value;
        
        title(ax_int, app.popup_int_spike.Items{int_idxs(i)})
        
        res_sm= nan(length(edges)-1,length(spks)-num_beats+2);
        peak = nan(1,length(spks)-num_beats+1);
        int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
        tmp = [];
        for j = 1: size(spks,2)
            tmp = [tmp,spks{1,j}]; 
        end

        [N,~] = histcounts(tmp,edges);
        res_sm(:,1) = movmean(N, mov_mean);

        for j = num_beats: size(spks,2)
            tmp = [];
            for k = -num_beats+1: 0
                tmp = [tmp,spks{1,j+k}];
            end
            [N,~] = histcounts(tmp,edges);
            res_sm(:,j-num_beats+2) = movmean(N, mov_mean);
            [~,idx] = max(res_sm(int(1):int(2),j-num_beats+1));
            peak(j-num_beats+1)= edges(int(1)+idx);
        
        end


        hb_ts = app.hb_res.t_events(app.hb_res.use_beats(:,int_idx));
        p=pcolor(ax_int,hb_ts(num_beats:end), edges(2:end),res_sm(:,2:end));
        p.EdgeColor = 'interp';
        line (ax_int,[hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
        line (ax_int,hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
        app.spike_res.spike_lag.line =line (app.ax_lag,[-1, -1], [edges(2), edges(end)], 'Color', 'r','LineStyle',':');
        xlim (ax_int,[hb_ts(num_beats),hb_ts(end)])

         col_idx =  rem(i,14);
        plts(end+1) = plot (ax_all, edges(2:end), res_sm(:,1), 'Color', cols{col_idx});
        xlim (ax_all, [edges(2),edges(end)])
        [~,idx_max] = max(res_sm(int(1):int(2),1));
        tmp = edges(int(1):int(2));
        lag = mean([tmp(idx_max),tmp(idx_max+1)]);
        yl = [min(res_sm(:,1)), max(res_sm(:,1))];
        line (ax_all, [lag,lag],yl , 'Color', cols{col_idx})
        leg{end+1} =  app.popup_int_spike.Items{int_idxs(i)};
        if i == length(int_idxs)
            legend(ax_all, plts,leg)
        end
        title (ax_all,'spike lag histogramms') 

        plts_n(end+1) = plot (ax_all_n, edges(2:end), res_sm(:,1)./std(res_sm(:,1)), 'Color', cols{col_idx});
        xlim (ax_all_n, [edges(2),edges(end)])
        yl = [min(res_sm(:,1)./std(res_sm(:,1))), max(res_sm(:,1)./std(res_sm(:,1)))];
        line (ax_all_n, [lag,lag],yl , 'Color', cols{col_idx})
        
        if i == length(int_idxs)
            legend(ax_all_n, plts_n, leg)
        end
        title (ax_all_n, 'normalized spike lag histogramms')

        title(ax_int,[leg{end} ' global lag: ' num2str(lag) ' [s]'  ])
        
        for k = 1 : length(form_idxs)
            switch form_idxs(k)
                case 1
                    h.Visible = 'on';
                    savefig(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_' app.popup_int_spike.Items{int_idx} '.fig'],'compact')
                    %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.fig'])
                case 2
                    saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_' app.popup_int_spike.Items{int_idx} '.jpeg'])
                case 3
                    saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_' app.popup_int_spike.Items{int_idx} '.epsc'])
            end
        end
        close(h)
    end

    
    for k = 1 : length(form_idxs)
            switch form_idxs(k)
                case 1
                    h_all.Visible = 'on';
                    savefig(h_all,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist.fig'],'compact')
                    %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.fig'])
                case 2
                    saveas(h_all,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist.jpeg'])
                case 3
                    saveas(h_all,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist.epsc'])
            end
    end
    close(h_all)
    
    for k = 1 : length(form_idxs)
            switch form_idxs(k)
                case 1
                    h_all_n.Visible = 'on';
                    savefig(h_all_n,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist_norm.fig'],'compact')
                    %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_rates_unit_' num2str(i) '.fig'])
                case 2
                    saveas(h_all_n,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist_norm.jpeg'])
                case 3
                    saveas(h_all_n,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_lag_hist_norm.epsc'])
            end
    end
    close(h_all_n)
end

printcell = cell([]);
for i = 1: max(app.spike_res.cluster)
    printcell{1,1,i} = ['Cluster ' num2str(i)];
    printcell{3,1,i} = 'Mean spike rate [Hz]';
    
    printcell{4,1,i} = 'empty cardiac cycles [%]';
    printcell{5,1,i} = 'PLV first';
    printcell{6,1,i} = 'PLV average';
    printcell{7,1,i} = 'PLV last';

    printcell{8,1,i} = 'phase linear first [rad]';
    printcell{9,1,i} = 'phase linear average [rad]';
    printcell{10,1,i} = 'phase linear last [rad]';

    printcell{11,1,i} = 'phase circular first [rad]';
    printcell{12,1,i} = 'phase circular average [rad]';
    printcell{13,1,i} = 'phase circular last [rad]';
    for j = 1: length(int_idxs)
        printcell{2,1+j,i} = simple_name(int_names{int_idxs(j)});
    end
    tmp = 1;
    printcell{15,tmp,i} = 'first spike phase lin [rad]';
    for j = 1: length(int_idxs)
        printcell{16,tmp,i} = simple_name(int_names{int_idxs(j)});
        tmp = tmp+1;
    end
    printcell{15,tmp,i} = 'mean spike phase lin [rad]';
    for j = 1: length(int_idxs)
        printcell{16,tmp,i} = simple_name(int_names{int_idxs(j)});
        tmp = tmp+1;
    end
    printcell{15,tmp,i} = 'last spike phase lin [rad]';
    for j = 1: length(int_idxs)
        printcell{16,tmp,i} = simple_name(int_names{int_idxs(j)});
        tmp = tmp+1;
    end
   
end
for i = 1: length(int_idxs)

    if int_idxs(i) ==1
        xl = ts;
        strt_stp = [1, size(data,1)];
    else
        xl = [min(min(app.burst_ints(int_idxs(i)-1).borders)), max(max(app.burst_ints(int_idxs(i)-1).borders))];
        strt_stp = [find(data(:,2) >= xl(1),1), find(data(:,2) >= xl(2),1)];
    end

    n_cluster = max(app.spike_res.cluster);
    cols = distinguishable_colors(n_cluster,[1 1 1; 0 0 0]);
    h = figure('Position', get(0, 'Screensize'),'Visible', 'off'); %%%%
    subplot (4,1,1:2)
    fs = 1/app.data(app.settings.channel_idx.msna).ts(1);           %%
    w0 = [300,5000]/(0.5*fs);
    if w0(2) >=1
        w0(2) = 0.999999;
    end
    [b,a] = butter(3,w0,"bandpass");               %%
    msna_int_f = filtfilt(b,a,data(strt_stp(1):strt_stp(2),1) );    %%
    msna_int_f(:,2) = data(strt_stp(1):strt_stp(2),2);

    leg_str = {[]};
    for j = 1: n_cluster        
        tmp = app.spike_res.spike_ts((app.spike_res.cluster == j) & app.spike_res.use_spikes(:,int_idxs(i))')/1000;
        plotspikes{j,1} = nan(length(tmp)*3,1);
        tmp_idx = 1:3:length(tmp)*3;
        plotspikes{j,1}(tmp_idx,1)= tmp;
        plotspikes{j,1}(tmp_idx+1,1)= tmp;
    
        plotspikes{j,2} = nan(length(tmp)*3,1);
        plotspikes{j,2}(tmp_idx,2) = -50;
        plotspikes{j,2}(tmp_idx+1,2)= 50;
        plotspikes{j,3} = nan(length(tmp)*3,1);
        plotspikes{j,3}(tmp_idx,2) = j -0.45;
        plotspikes{j,3}(tmp_idx+1,2) = j +0.45;
        
        subplot (4,1,1:2)
        hold on 
        plot (plotspikes{j,1},plotspikes{j,3},'Color', cols(j,:), 'LineWidth', 0.7 )
        hold off
        xlabel('[s]')
        ylabel ('Cluster')
        subplot (4,1,4)
        hold on 
        plot ((-20:1:20)*app.data(app.settings.channel_idx.msna).ts(1),mean(app.spike_res.spikes(:,app.spike_res.cluster ==j),2),'Color', cols(j,:), 'LineWidth',1.5 )
        leg_str{j} = ['n = ' num2str(sum(app.spike_res.cluster ==j))];
        hold off
        xlabel('[s]')

    end

    subplot (4,1,1:2)
    xlim (xl)
    ylim([.5 n_cluster+.5])
    title(['Spike timestamps ' simple_name(int_names{int_idxs(i)}) ' all units'])
    
    subplot (4,1,4)
    legend(leg_str,'Location','southwest','FontSize',7)
    title(['Average spikeshape ' simple_name(int_names{int_idxs(i)}) ' all units'])
    
    subplot(4,1,3)
    if ~isempty(app.hb_res)
        hb = nan(length(app.hb_res.t_events(app.hb_res.use_beats(:,int_idxs(i))))*3,2);
        tmp_idx = 1:3:length(app.hb_res.t_events(app.hb_res.use_beats(:,int_idxs(i))))*3;
        hb(tmp_idx,1)= app.hb_res.t_events(app.hb_res.use_beats(:,int_idxs(i)));
        hb(tmp_idx+1,1)= app.hb_res.t_events(app.hb_res.use_beats(:,int_idxs(i)));
        hb(tmp_idx,2)= -50;
        hb(tmp_idx+1,2)= 50;
        yyaxis left
        hold off
        plot (hb(:,1),hb(:,2),'LineWidth',1.5,'Color','k')
        hold on
    end
    plot( msna_int_f(:,2),  msna_int_f(:,1), 'Color',[0.6 0.7 0.7],'LineStyle','-')
    %plot(downsample((1:length(data))*app.data(app.settings.channel_idx.msna).ts(1),100),downsample(data(:,1),100), 'Color',[0.2 0.2 0.3],'LineStyle','-') %%%%% only plot msna in limits
    for j = 1: n_cluster 
        plot(plotspikes{j,1},plotspikes{j,2},'Color', cols(j,:), 'LineWidth', 0.5 )
    end
    ylim([-20 20])
    xlim(xl)
    hold off
    title(['MSNA ' int_names{int_idxs(i)} ' all units'])
    %% better spikefreq
    spikemap = zeros(length(data),1);
    spikemap(app.spike_res.spike_idx) = 1;
    m = downsample(movsum(spikemap, 5000),20);
    
    %%
    
    yyaxis right 
    hold off
    plot((1:length(m))*app.data(app.settings.channel_idx.msna).ts(1)*20,m,'Color',[0 0.7 0],'LineWidth',1.4,'Marker','none');
    xlim(xl)
    ylim([min(m) max(m)])

    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                h.Visible = 'on';
                savefig(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_overview_' simple_name(int_names{int_idxs(i)}) '.fig'],'compact')
                %switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_overview_' simple_name(int_names{int_idxs(i)}) '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_overview_' simple_name(int_names{int_idxs(i)}) '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_overview_' simple_name(int_names{int_idxs(i)}) '.epsc'])
        end
    end
    close(h)
    ndig = length(num2str(n_cluster));
    num_templ = [];
    for j = 1 :ndig
        num_templ = [num_templ '0'];
    end

    for j = 1: n_cluster
        num_str = num_templ;
        tmp_num = num2str(j);
        num_str(end-length(tmp_num)+1:end) = tmp_num;
        h2 = figure('Position', get(0, 'Screensize'),'Visible', 'off');
        subplot (4,3,1:3)
        
        yyaxis left
       
        hold off 
        if ~isempty(app.hb_res)
            plot (hb(:,1),hb(:,2),'LineWidth',1.5,'Color','k')
        end
        hold on
        plot( msna_int_f(:,2),  msna_int_f(:,1), 'Color',[0.6 0.7 0.7],'LineStyle','-')
        %plot(downsample((1:length(data))*app.data(app.settings.channel_idx.msna).ts(1),100),downsample(data(:,1),100), 'Color',[0.2 0.2 0.3],'LineStyle','-') %%% only plot msna in interval
        plot(plotspikes{j,1},plotspikes{j,2},'Color', cols(j,:), 'LineWidth', 0.5)
        ylim([-20 20])
        xlim(xl)
        hold off
        title(['MSNA ' simple_name(int_names{int_idxs(i)}) ' unit: ' num2str(j)])
        yyaxis right 
        hold off
        
spikemap = zeros(length(data),1);                                   %% check
spikemap(app.spike_res.spike_idx(app.spike_res.cluster == j)) = 1;  %% check
m = downsample(movsum(spikemap, 5000),20);                          %% check
        
        plot((1:length(m))*app.data(app.settings.channel_idx.msna).ts(1)*20,m,'LineWidth',1.4,'Marker','none');%,'Color',[0 0.7 0]
        xlim(xl)
        ylim([min(m) max(m)])
        %% only with HB
        if ~isempty(app.hb_res)
            subplot (4,3,4)
            
            tmp_idx = app.spike_res.spike(int_idxs(i)).first{j}(:,2);
            tmp = mean(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),2);
            tmp_s = std(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),0,2);
            
            hold off
            fill([(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'r','FaceAlpha',0.3,'EdgeColor','none')
            hold on
            plot((-20:1:20) /10, tmp, 'r', 'LineWidth',2)
            hold off
            title(['first spikes Unit: ' num2str(j)])
            
            subplot (4,3,5)
    
            tmp_idx = app.spike_res.spike(int_idxs(i)).last{j}(:,2);
            tmp = mean(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),2);
            tmp_s = std(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),0,2);
            
            hold off
            fill([(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'b','FaceAlpha',0.3,'EdgeColor','none')
            hold on
            plot((-20:1:20) /10, tmp, 'b', 'LineWidth',2)
            hold off
            title(['all spikes Unit: ' num2str(j)])
            
            subplot (4,3,6)
            
            tmp = mean(app.spike_res.spikes(:,app.spike_res.cluster == j),2);
            tmp_s = std(app.spike_res.spikes(:,app.spike_res.cluster == j),0,2);
            
            hold off
            fill([(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'g','FaceAlpha',0.3,'EdgeColor','none')
            hold on
            plot((-20:1:20) /10, tmp, 'g', 'LineWidth',2)
            hold off
            title(['last spikes Unit: ' num2str(j)])
    
            subplot (4,3,7)
    
            [aaa, r, phi, zmzm] = circ_plot_GDA(app.spike_res.spike(int_idxs(i)).totXrad_min{j}','pretty','ro',true,'linewidth',4,'Color','r');
            if phi<0
                phi=2*pi-abs(phi);
            end
            
            circFIRST=mean(phi);
            plvFIRST=mean(r);
            
            hold on
            plot([-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
            plot([-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
            plot([-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
            plot([0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
            dummyh = line(nan, nan, 'Marker', 'none', 'Color', 'r');
            legend(dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
            printcell{5,i+1,j} = plvFIRST;
            printcell{3,i+1,j} = app.spike_res.spike(int_idxs(i)).m_fr(j);
            printcell{4,i+1,j} = length(find(isnan(app.spike_res.spike(int_idxs(i)).first{1, j}(:,1))))/size(app.spike_res.spike(int_idxs(i)).first{1, j}(:,1),1)*100;

            tmp =  num2cell(app.spike_res.spike(int_idxs(i)).totXrad_min{1, j}');

            printcell(17:16+length(tmp),i,j) = tmp;
            title('FIRST Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_min{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
            printcell{8,i+1,j} = mean(app.spike_res.spike(int_idxs(i)).totXrad_min{j});
            printcell{11,i+1,j} = mean(phi);
            set(gca,'XColor', 'none','YColor','none')
            hold off
    
            subplot (4,3,8)
    
            [aaa, r, phi, zmzm] = circ_plot_GDA(app.spike_res.spike(int_idxs(i)).totXrad_mean{j}','pretty','go',true,'linewidth',4,'Color','g');
            if phi<0
                phi=2*pi-abs(phi);
            end
            
            circFIRST=mean(phi);
            plvFIRST=mean(r);
            
            hold on
            plot([-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
            plot([-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
            plot([-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
            plot([0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
            dummyh = line(nan, nan, 'Marker', 'none', 'Color', 'g');
            legend(dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
            printcell{6,i+1,j} = plvFIRST;
            tmp =  num2cell(app.spike_res.spike(int_idxs(i)).totXrad_mean{1, j}');
            printcell(17:16+length(tmp),i+length(int_idxs),j) = tmp;
            title('AVERAGE Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_mean{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
            printcell{9,i+1,j} = mean(app.spike_res.spike(int_idxs(i)).totXrad_mean{j});
            printcell{12,i+1,j} = mean(phi);
            set(gca, 'XColor', 'none','YColor','none')
            hold off
    
            subplot (4,3,9)
    
            [aaa, r, phi, zmzm] = circ_plot_GDA(app.spike_res.spike(int_idxs(i)).totXrad_max{j}','pretty','bo',true,'linewidth',4,'Color','b');
            if phi<0
                phi=2*pi-abs(phi);
            end
            
            circFIRST=mean(phi);
            plvFIRST=mean(r);
            
            hold on
            plot([-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
            plot([-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
            plot([-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
            plot([0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
            dummyh = line(nan, nan, 'Marker', 'none', 'Color', 'b');
            legend(dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
            printcell{7,i+1,j} = plvFIRST;
            
            tmp =  num2cell(app.spike_res.spike(int_idxs(i)).totXrad_max{1, j}');
            printcell(17:16+length(tmp),i+2*length(int_idxs),j) = tmp;
            title('LAST Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_max{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
            printcell{10,i+1,j} = mean(app.spike_res.spike(int_idxs(i)).totXrad_max{j});
            printcell{13,i+1,j} = mean(phi);
            set(gca,'XColor', 'none','YColor','none')
            hold off
    
            subplot (4,3,10)
    
            nhist(app.spike_res.spike(int_idxs(i)).totXrad_min{j},'color',[.8 .3 .3],'text','box','median','noerror');
            xlabel('First Spike Phase [rad]') 
            hold on
            xlim([0 2*pi])
            hold off
            grid on 
            grid minor
    
            subplot (4,3,11)
    
            nhist(app.spike_res.spike(int_idxs(i)).totXrad_mean{j},'color',[.3 .8 .3],'text','box','median','noerror');
            xlabel('Average Firing Phase [rad]') 
            hold on
            xlim([0 2*pi])
            hold off
            grid on
            grid minor
    
            subplot (4,3,12)
    
            nhist(app.spike_res.spike(int_idxs(i)).totXrad_max{j},'color',[.3 .3 .8],'text','box','median','noerror');
            xlabel('Last Spike Phase [rad]') 
            hold on
            xlim([0 2*pi])
            hold off
            grid on
            grid minor
        end
        for k = 1 : length(form_idxs)
            switch form_idxs(k)
                case 1
                    h2.Visible = 'on';
                    savefig(h2,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_unit_' num_str '_' simple_name(int_names{int_idxs(i)}) '.fig'],'compact')
                    % switch_vis([path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_unit_' num_str '_' simple_name(int_names{int_idxs(i)}) '.fig'])
                case 2
                    saveas(h2,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_unit_' num_str '_' simple_name(int_names{int_idxs(i)}) '.jpeg'])
                case 3
                    saveas(h2,[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_spike_analysis_unit_' num_str '_' simple_name(int_names{int_idxs(i)}) '.epsc'])
            end
        end
    close(h2)
    app.lbl_working.Text = [num2str(round( (100/length(int_idxs))*(i-1) + ((100/length(int_idxs))/n_cluster)*j )) '% done'];
    drawnow %pause(.01) %%%%
    end
end
for i =1: size(printcell,3)
    writecell(printcell(:,:,i),[path '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) 'spike_results.xls'], 'Sheet', printcell{1,1,i})
    
end
end