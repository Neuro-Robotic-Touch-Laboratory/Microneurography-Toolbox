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
    bin_step = app.edt_step.Value;
    nbins = app.edt_max.Value;
    cols = distinguishable_colors(max(app.spike_res.cluster),[1 1 1; 1 0 0; 0 0 1; 0.80 0.95 0.90]);
    for i= 1:max(app.spike_res.cluster)
        
        h = figure('Position', get(0, 'Screensize'),'Visible', 'off');
        
        subplot(1,2,1)
        mrate = nan(length(app.lstbx_ints.Value),1);
        
        for j = 1 :length(app.lstbx_ints.Value)
            tmp_idx = find(strcmp(app.lstbx_ints.Value{j},app.lstbx_ints.Items));
            mrate(j) = app.spike_res.spike(tmp_idx).m_fr(i);
        end
        ax = gca;
        
        plot (ax,mrate, 'LineWidth', 2, 'LineStyle','--', 'Color',cols(i,:),'Marker','diamond','MarkerSize',10,'MarkerFaceColor',cols(i,:))
        
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
            title('FIRST Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_min{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
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
            title('FIRST Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_min{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
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
            title('FIRST Spike Phase:' );
            subtitle(['linear=' num2str(mean(app.spike_res.spike(int_idxs(i)).totXrad_min{j})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
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

end