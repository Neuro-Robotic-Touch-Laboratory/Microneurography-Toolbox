function update_spike_axis(app, static)
%UPDATE_SPIKE_AXIS Summary of this function goes here
%   Detailed explanation goes here


sel_cluster = find(strcmp(app.popup_cluster.Value,app.popup_cluster.Items));
int_idx = find(strcmp(app.popup_int_spike.Value,app.popup_int_spike.Items));

%data = app.data(app.settings.channel_idx.msna).data(app.settings.interval(1,1)/app.data(app.settings.channel_idx.msna).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.msna).ts(1));
[data,ts,~, ~] = current_signal(app, app.settings.channel_idx.msna);

if int_idx ==1
    xl = ts;
else
    xl = [min(min(app.burst_ints(int_idx-1).borders)), max(max(app.burst_ints(int_idx-1).borders))];
end

n_cluster = max(app.spike_res.cluster);
cols = distinguishable_colors(n_cluster,[1 1 1; 0 0 0]);
 %% plot all spike of one cluster as one object
plotspikes = {[]}; 
if static
    cla(app.ax_cluster_timestamps)
    hold(app.ax_cluster_timestamps, 'on')
    cla(app.ax_cluster_shapes)
    hold(app.ax_cluster_shapes, 'on')
end
leg_str = {[]};
for i = 1: n_cluster        
    tmp = app.spike_res.spike_ts((app.spike_res.cluster == i) & app.spike_res.use_spikes(:,int_idx)')/1000;
    plotspikes{i,1} = nan(length(tmp)*3,1);
    tmp_idx = 1:3:length(tmp)*3;
    plotspikes{i,1}(tmp_idx,1)= tmp;
    plotspikes{i,1}(tmp_idx+1,1)= tmp;

    plotspikes{i,2} = nan(length(tmp)*3,1);
    plotspikes{i,2}(tmp_idx,2) = -50;
    plotspikes{i,2}(tmp_idx+1,2)= 50;
    plotspikes{i,3} = nan(length(tmp)*3,1);
    plotspikes{i,3}(tmp_idx,2) = i -0.45;
    plotspikes{i,3}(tmp_idx+1,2) = i +0.45;

    if static
        plot (app.ax_cluster_timestamps,plotspikes{i,1},plotspikes{i,3},'Color', cols(i,:), 'LineWidth', 0.7 )
        plot (app.ax_cluster_shapes,(-20:1:20)*app.data(app.settings.channel_idx.msna).ts(1),app.spike_res.spike(int_idx).shape{i},'Color', cols(i,:), 'LineWidth',1.5 )
        leg_str{i} = ['n = ' num2str(app.spike_res.spike(int_idx).n_spikes(i))];
        xlim (app.ax_cluster_timestamps,xl)
        
    end
end

if static
    lgnd = legend(app.ax_cluster_shapes,leg_str,'Location','southwest','FontSize',7);
    lgnd.ItemTokenSize(1) = 8;
    hold(app.ax_cluster_timestamps, 'off')
    hold(app.ax_cluster_shapes, 'off')
end
hold(app.ax_cluster_timestamps, 'off')
ylim(app.ax_cluster_timestamps,[.5 n_cluster+.5])

%% plot all heartbeats as one object

hb = nan(length(app.hb_res.t_events(app.hb_res.use_beats(:,int_idx)))*3,2);
tmp_idx = 1:3:length(app.hb_res.t_events(app.hb_res.use_beats(:,int_idx)))*3;
hb(tmp_idx,1)= app.hb_res.t_events(app.hb_res.use_beats(:,int_idx));
hb(tmp_idx+1,1)= app.hb_res.t_events(app.hb_res.use_beats(:,int_idx));
hb(tmp_idx,2)= -50;
hb(tmp_idx+1,2)= 50;
yyaxis(app.ax_spike_msna, 'left')
cla(app.ax_spike_msna)
hold(app.ax_spike_msna, 'off')
plot (app.ax_spike_msna,hb(:,1),hb(:,2),'LineWidth',1.5,'Color','k')
hold(app.ax_spike_msna, 'on')

plot(app.ax_spike_msna, downsample((1:length(data))*app.data(app.settings.channel_idx.msna).ts(1),100),downsample(data,100),'LineStyle','-')%, 'Color',[0.2 0.2 0.3]
plot (app.ax_spike_msna,plotspikes{sel_cluster,1},plotspikes{sel_cluster,2},'Color', cols(sel_cluster,:), 'LineWidth', 0.7 ,'LineStyle','-.' )
ylim(app.ax_spike_msna,[-20 20])
xlim(app.ax_spike_msna,[0 length(data)*app.data(app.settings.channel_idx.msna).ts(1)])
hold(app.ax_spike_msna, 'off')

app.lbl_spike_ttl2.Text = ['Timestamps unit: ' num2str(sel_cluster)];
app.lbl_spike_ttl2.FontColor = cols(sel_cluster,:);
app.lbl_spike_ttl3.Text = ['Firing frequency unit: ' num2str(sel_cluster)];

%% better spikefreq
spikemap = zeros(length(data),1);
spikemap(app.spike_res.spike_idx(app.spike_res.cluster == sel_cluster)) = 1;
m = downsample(movsum(spikemap, 5000),20);

%%

yyaxis(app.ax_spike_msna, 'right') 
hold(app.ax_spike_msna, 'off')
cla(app.ax_spike_msna)
plot(app.ax_spike_msna,(1:length(m))*app.data(app.settings.channel_idx.msna).ts(1)*20,m,'LineWidth',1.4,'Marker','none');%,'Color',[0 0.7 0]
xlim(app.ax_spike_msna, xl)
app.spike_res.rate_ylim = [min(m) max(m)];
ylim(app.ax_spike_msna,app.spike_res.rate_ylim)
yyaxis(app.ax_spike_msna, 'left')
ylim(app.ax_spike_msna,[-20 20])


tmp_idx = app.spike_res.spike(int_idx).first{sel_cluster}(:,2);
tmp = mean(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),2);
tmp_s = std(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),0,2);

hold(app.ax_1st_shape, 'off')
fill(app.ax_1st_shape,[(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'r','FaceAlpha',0.3,'EdgeColor','none')
hold(app.ax_1st_shape, 'on')
plot(app.ax_1st_shape, (-20:1:20) /10, tmp, 'r', 'LineWidth',2)
hold(app.ax_1st_shape, 'off')

tmp_idx = app.spike_res.spike(int_idx).last{sel_cluster}(:,2);
tmp = mean(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),2);
tmp_s = std(app.spike_res.spikes(:,tmp_idx(~isnan(tmp_idx))),0,2);

hold(app.ax_last_shape, 'off')
fill(app.ax_last_shape,[(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'b','FaceAlpha',0.3,'EdgeColor','none')
hold(app.ax_last_shape, 'on')
plot(app.ax_last_shape, (-20:1:20) /10, tmp, 'b', 'LineWidth',2)
hold(app.ax_last_shape, 'off')

tmp = mean(app.spike_res.spikes(:,app.spike_res.cluster == sel_cluster),2);
tmp_s = std(app.spike_res.spikes(:,app.spike_res.cluster == sel_cluster),0,2);

hold(app.ax_all_shape, 'off')
fill(app.ax_all_shape,[(-20:1:20),(20:-1:-20)] /10, [tmp'+tmp_s' , flip(tmp'-tmp_s')], 'g','FaceAlpha',0.3,'EdgeColor','none')
hold(app.ax_all_shape, 'on')
plot(app.ax_all_shape, (-20:1:20) /10, tmp, 'g', 'LineWidth',2)
hold(app.ax_all_shape, 'off')



[aaa, r, phi, zmzm]=circ_plot_ax(app.ax_1st_phase, app.spike_res.spike(int_idx).totXrad_min{sel_cluster}','pretty','ro',true,'linewidth',4,'Color','r');
if phi<0
    phi=2*pi-abs(phi);
end

circFIRST=mean(phi);
plvFIRST=mean(r);

hold(app.ax_1st_phase,'on')
plot(app.ax_1st_phase,[-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
plot(app.ax_1st_phase,[-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
plot(app.ax_1st_phase,[-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
plot(app.ax_1st_phase,[0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
dummyh = line(app.ax_1st_phase,nan, nan, 'Marker', 'none', 'Color', 'r');
legend(app.ax_1st_phase,dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
%title(app.ax_1st_phase,'FIRST Spike Phase:' );
app.lbl_1st_phase_lin.Text = ['linear = ' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad'];
app.lbl_1st_phase_circ.Text = [ 'circular = ' num2str(mean(phi)) ' rad'];
%subtitle(app.ax_1st_phase,['linear=' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
set(app.ax_1st_phase,'XColor', 'none','YColor','none')
hold(app.ax_1st_phase,'off')


[aaa, r, phi, zmzm] = circ_plot_ax(app.ax_all_phase, app.spike_res.spike(int_idx).totXrad_mean{sel_cluster}','pretty','go',true,'linewidth',4,'Color','g');
if phi<0
    phi=2*pi-abs(phi);
end

circFIRST=mean(phi);
plvFIRST=mean(r);

hold(app.ax_all_phase,'on')
plot(app.ax_all_phase,[-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
plot(app.ax_all_phase,[-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
plot(app.ax_all_phase,[-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
plot(app.ax_all_phase,[0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
dummyh = line(app.ax_all_phase,nan, nan, 'Marker', 'none', 'Color', 'g');
legend(app.ax_all_phase,dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
% title(app.ax_all_phase,'FIRST Spike Phase:' );
app.lbl_all_phase_lin.Text = ['linear = ' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad'];
app.lbl_all_phase_circ.Text = [ 'circular = ' num2str(mean(phi)) ' rad'];
% subtitle(app.ax_all_phase,['linear=' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
set(app.ax_all_phase,'XColor', 'none','YColor','none')
hold(app.ax_all_phase,'off')

[aaa, r, phi, zmzm] = circ_plot_ax(app.ax_last_phase, app.spike_res.spike(int_idx).totXrad_max{sel_cluster}','pretty','bo',true,'linewidth',4,'Color','b');
if phi<0
    phi=2*pi-abs(phi);
end

circFIRST=mean(phi);
plvFIRST=mean(r);

hold(app.ax_last_phase,'on')
plot(app.ax_last_phase,[-sqrt(2)/2 sqrt(2)/2],[sqrt(2)/2 -sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % plot
plot(app.ax_last_phase,[-sqrt(2)/2 sqrt(2)/2],[-sqrt(2)/2 sqrt(2)/2],'--','Color',[0.2 0.2 0.2])     % grid
plot(app.ax_last_phase,[-1 1],[0 0],'--','Color',[0.2 0.2 0.2])                                      % lines
plot(app.ax_last_phase,[0 0],[-1 1],'--','Color',[0.2 0.2 0.2])                                      % 
dummyh = line(app.ax_last_phase,nan, nan, 'Marker', 'none', 'Color', 'b');
legend(app.ax_last_phase,dummyh, ['PLV= ' num2str(mean(r))], 'Location','best')
% title(app.ax_last_phase,'FIRST Spike Phase:' );
app.lbl_last_phase_lin.Text = ['linear = ' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad'];
app.lbl_last_phase_circ.Text = [ 'circular = ' num2str(mean(phi)) ' rad'];
% subtitle(app.ax_last_phase,['linear=' num2str(mean(app.spike_res.spike(int_idx).totXrad_min{sel_cluster})) ' rad' '           circular= ' num2str(mean(phi)) ' rad'])
set(app.ax_last_phase,'XColor', 'none','YColor','none')
hold(app.ax_last_phase,'off')

%% plot hist 1st spikes
nhist_ax(app.spike_res.spike(int_idx).totXrad_min{sel_cluster},'color',[.8 .3 .3],'text','box','median','noerror','axis', app.ax_1st_hist)
xlabel(app.ax_1st_hist,'First Spike Phase [rad]') 
hold(app.ax_1st_hist, 'on')
xlim(app.ax_1st_hist,[0 2*pi])
hold(app.ax_1st_hist, 'off')
grid(app.ax_1st_hist, 'on') 
grid(app.ax_1st_hist, 'minor')


%% plott his all spikes
nhist_ax(app.spike_res.spike(int_idx).totXrad_mean{sel_cluster},'color',[.3 .8 .3],'text','box','median','noerror','axis', app.ax_all_hist)
xlabel(app.ax_all_hist,'Average Firing Phase [rad]') 
hold(app.ax_all_hist, 'on')
xlim(app.ax_all_hist, [0 2*pi])
hold(app.ax_all_hist, 'off')
grid(app.ax_all_hist, 'on')
grid(app.ax_all_hist, 'minor')


%% plot hist last spikes
nhist_ax(app.spike_res.spike(int_idx).totXrad_max{sel_cluster},'color',[.3 .3 .8],'text','box','median','noerror','axis', app.ax_last_hist)
xlabel(app.ax_last_hist, 'Last Spike Phase [rad]') 
hold(app.ax_last_hist, 'on')
xlim(app.ax_last_hist,[0 2*pi])
hold(app.ax_last_hist, 'off')
grid(app.ax_last_hist, 'on') 
grid(app.ax_last_hist, 'minor')

ylim(app.ax_cluster_timestamps,[.5 n_cluster+.5])
end

