function update_burst_axis(app, keep_xl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cla(app.ax_msna_int), cla(app.ax_burst_res)  
cla(app.ax_burst_interval_his), cla(app.ax_burst_duration_his)
cla(app.ax_burst_amplitude_his), cla(app.ax_burst_integral_his)
cla(app.ax_burst_res), cla(app.ax_msna_int)

t_msna = app.burst_res.ts(1) : app.burst_res.ts(2) : app.burst_res.ts(3);
t_i = app.settings.interval(1,1);
t_f = app.settings.interval(1,2);
int_idx = find(strcmp(app.popup_int.Value, app.popup_int.Items));
if keep_xl
    xl = xlim(app.ax_msna_int);
else
    if int_idx ==1
    xl = [app.burst_res.ts(1), app.burst_res.ts(3)];
    else
        xl = [min(min(app.burst_ints(int_idx-1).borders)), max(max(app.burst_ints(int_idx-1).borders))];
    end
end
yl = ylim(app.ax_msna_raw);

hold(app.ax_msna_raw,'off')
for i = 1: size(app.settings.burst_rem_int,1)
    tmp = app.settings.burst_rem_int(i,:);
    fill(app.ax_msna_raw,[tmp(1), tmp(1), tmp(2), tmp(2)],[yl(1), yl(2), yl(2), yl(1)],'r','FaceAlpha', 0.3, 'HitTest','off')
    hold(app.ax_msna_raw,'on')
end

plot(app.ax_msna_raw, t_msna,app.data(int32(app.settings.channel_idx.msna)).data(int32(t_i/app.burst_res.ts(2)+1):int32(t_f/app.burst_res.ts(2))),'LineWidth',1.2)
yl = quantile(app.data(int32(app.settings.channel_idx.msna)).data(int32(t_i/app.burst_res.ts(2)+1):int32(t_f/app.burst_res.ts(2))),[0.00001 .99999]);
ylim(app.ax_msna_raw,yl)
xlim(app.ax_msna_raw,[t_msna(1) t_msna(end)])

if ~isempty(app.hb_res)
    hb_idx = app.hb_res.use_beats(:,1) & app.hb_res.use_beats(:,int_idx);
end
burst_idx = app.burst_res.use_burst(:,1) & app.burst_res.use_burst(:,2) & app.burst_res.use_burst(:,int_idx+1);
if int_idx ==1
    dur = diff([app.burst_res.ts(1),app.burst_res.ts(end)]);
else
    dur = sum(diff(app.burst_ints(int_idx-1).borders,1,2));
end

yl = [min(app.burst_res.x(2001:end)) max(app.burst_res.x(2001:end))];
hold(app.ax_msna_int,'off')
for i = 1: size(app.settings.burst_rem_int,1)
    tmp = app.settings.burst_rem_int(i,:);
    fill(app.ax_msna_int,[tmp(1), tmp(1), tmp(2), tmp(2)],[yl(1), yl(2), yl(2), yl(1)],'r','FaceAlpha', 0.3, 'HitTest','off')
    hold(app.ax_msna_int,'on')
end

plot(app.ax_msna_int, t_msna(2001:end), app.burst_res.x(2001:end),'LineWidth',1.2,'HitTest','off')  
ylim(app.ax_msna_int, yl)  
hold(app.ax_msna_int, 'on')


%% plot all heartbeats as one object 
if ~isempty(app.hb_res)
    hb = nan(length(app.hb_res.t_events(hb_idx))*3,2);
    tmp_idx = 1:3:length(app.hb_res.t_events(hb_idx))*3;
    hb(tmp_idx,1)= app.hb_res.t_events(hb_idx);
    hb(tmp_idx+1,1)= app.hb_res.t_events(hb_idx);
    yl = ylim(app.ax_msna_int);
    hb(tmp_idx,2)= yl(1);
    hb(tmp_idx+1,2)= yl(2);
    plot (app.ax_msna_int,hb(:,1),hb(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2],'HitTest','off')
end
xlim(app.ax_msna_int,xl)

%%
yl = ylim(app.ax_burst_res);
hold(app.ax_burst_res,'off')
for i = 1: size(app.settings.burst_rem_int,1)
    tmp = app.settings.burst_rem_int(i,:);
    fill(app.ax_burst_res,[tmp(1), tmp(1), tmp(2), tmp(2)],[yl(1), yl(2), yl(2), yl(1)],'r','FaceAlpha', 0.3, 'HitTest','off')
    hold(app.ax_burst_res,'on')
end

hold(app.ax_msna_int, 'on')
tmp=find(burst_idx);
burst_rem_idx = ~app.burst_res.use_burst(:,1) & app.burst_res.use_burst(:,2);
tmp_rem = find(burst_rem_idx);
        
colors = [1.0000, 0, 0;...
          0, 1.0000, 0;...
          0, 0, 0.1724;...
          1.0000, 0.1034, 0.7241;...
          1.0000, 0.8276, 0];

burst_plot1 =[nan,nan]; burst_plot2 =[nan,nan]; burst_plot3 =[nan,nan]; burst_plot4 =[nan,nan]; burst_plot5 =[nan,nan]; burst_plot_rem = [nan,nan];
burst_amp_plot1 = [nan,nan]; burst_amp_plot2 = [nan,nan]; burst_amp_plot3 = [nan,nan]; burst_amp_plot4 = [nan,nan]; burst_amp_plot5 = [nan,nan]; burst_amp_plot_rem = [nan,nan]; 

for j = 1:5: ceil(length(tmp)/5)*5 %%%% floor
    if j <= length(tmp)
        burst_plot1 = [burst_plot1;[(t_msna(app.burst_res.burst_loc(tmp(j),1):app.burst_res.burst_loc(tmp(j),2)))',...
                       (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp(j),1):app.burst_res.burst_loc(tmp(j),2)))'];...
                       [t_msna(app.burst_res.burst_loc(tmp(j),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp(j),1))];...
                       [nan,nan]];
        burst_amp_plot1 = [burst_amp_plot1 ; [app.burst_res.t_burst(tmp(j)) ,app.burst_res.burst_int(tmp(j))]];
    end
    if j+1 <= length(tmp)   
        burst_plot2 = [burst_plot2;[(t_msna(app.burst_res.burst_loc(tmp(j+1),1):app.burst_res.burst_loc(tmp(j+1),2)))',...
                       (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+1),1):app.burst_res.burst_loc(tmp(j+1),2)))'];...
                       [t_msna(app.burst_res.burst_loc(tmp(j+1),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+1),1))];...
                       [nan,nan]];
        burst_amp_plot2 = [burst_amp_plot2 ; [app.burst_res.t_burst(tmp(j+1)) ,app.burst_res.burst_int(tmp(j+1))]];
    end
    if j+2 <= length(tmp)
        burst_plot3 = [burst_plot3;[(t_msna(app.burst_res.burst_loc(tmp(j+2),1):app.burst_res.burst_loc(tmp(j+2),2)))',...
                       (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+2),1):app.burst_res.burst_loc(tmp(j+2),2)))'];...
                       [t_msna(app.burst_res.burst_loc(tmp(j+2),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+2),1))];...
                       [nan,nan]];
        burst_amp_plot3 = [burst_amp_plot3 ; [app.burst_res.t_burst(tmp(j+2)) ,app.burst_res.burst_int(tmp(j+2))]];
    end
    if j+3 <= length(tmp)
        burst_plot4 = [burst_plot4;[(t_msna(app.burst_res.burst_loc(tmp(j+3),1):app.burst_res.burst_loc(tmp(j+3),2)))',...
                       (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+3),1):app.burst_res.burst_loc(tmp(j+3),2)))'];...
                       [t_msna(app.burst_res.burst_loc(tmp(j+3),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+3),1))];...
                       [nan,nan]];
        burst_amp_plot4 = [burst_amp_plot4 ; [app.burst_res.t_burst(tmp(j+3)) ,app.burst_res.burst_int(tmp(j+3))]];
    end
    if j+4 <= length(tmp) 
        burst_plot5 = [burst_plot5;[(t_msna(app.burst_res.burst_loc(tmp(j+4),1):app.burst_res.burst_loc(tmp(j+4),2)))',...
                       (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+4),1):app.burst_res.burst_loc(tmp(j+4),2)))'];...
                       [t_msna(app.burst_res.burst_loc(tmp(j+4),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp(j+4),1))];...
                       [nan,nan]];
        burst_amp_plot5 = [burst_amp_plot5 ; [app.burst_res.t_burst(tmp(j+4)) ,app.burst_res.burst_int(tmp(j+4))]];
    end
end

for j = 1:length(tmp_rem) 
    burst_plot_rem = [burst_plot_rem;[(t_msna(app.burst_res.burst_loc(tmp_rem(j),1):app.burst_res.burst_loc(tmp_rem(j),2)))',...
                   (app.burst_res.htresh +app.burst_res.xshift(app.burst_res.burst_loc(tmp_rem(j),1):app.burst_res.burst_loc(tmp_rem(j),2)))'];...
                   [t_msna(app.burst_res.burst_loc(tmp_rem(j),1)), app.burst_res.htresh+app.burst_res.xshift(app.burst_res.burst_loc(tmp_rem(j),1))];...
                   [nan,nan]];
    burst_amp_plot_rem = [burst_amp_plot_rem ; [app.burst_res.t_burst(tmp_rem(j)) ,app.burst_res.burst_int(tmp_rem(j))]];
end

plot(app.ax_msna_int,burst_plot1(:,1),burst_plot1(:,2),'Color',colors(1,:),'LineWidth',1.5)
plot(app.ax_msna_int,burst_plot2(:,1),burst_plot2(:,2),'Color',colors(2,:),'LineWidth',1.5)
plot(app.ax_msna_int,burst_plot3(:,1),burst_plot3(:,2),'Color',colors(3,:),'LineWidth',1.5)
plot(app.ax_msna_int,burst_plot4(:,1),burst_plot4(:,2),'Color',colors(4,:),'LineWidth',1.5)
plot(app.ax_msna_int,burst_plot5(:,1),burst_plot5(:,2),'Color',colors(5,:),'LineWidth',1.5)
plot(app.ax_msna_int,burst_plot_rem(:,1),burst_plot_rem(:,2),'Color',[.7,.7,.7],'LineWidth',1.5)

hold(app.ax_burst_res,'on')
plot(app.ax_burst_res,burst_amp_plot1(:,1),burst_amp_plot1(:,2),'.','MarkerSize',20,'Color',colors(1,:))
plot(app.ax_burst_res,burst_amp_plot2(:,1),burst_amp_plot2(:,2),'.','MarkerSize',20,'Color',colors(2,:))
plot(app.ax_burst_res,burst_amp_plot3(:,1),burst_amp_plot3(:,2),'.','MarkerSize',20,'Color',colors(3,:))
plot(app.ax_burst_res,burst_amp_plot4(:,1),burst_amp_plot4(:,2),'.','MarkerSize',20,'Color',colors(4,:))
plot(app.ax_burst_res,burst_amp_plot5(:,1),burst_amp_plot5(:,2),'.','MarkerSize',20,'Color',colors(5,:))
plot(app.ax_burst_res,burst_amp_plot_rem(:,1),burst_amp_plot_rem(:,2),'.','MarkerSize',20,'Color',[.7,.7,.7])

plot(app.ax_burst_res, app.burst_res.t_burst(burst_idx), app.burst_res.burst_int(burst_idx), '-.','linewidth',2,'color','black','HitTest','off')
hold(app.ax_burst_res,'off')
 
hold(app.ax_burst_res, 'off')

xlim(app.ax_burst_res, xl)

CC.burst_integral = app.burst_res.burst_int(burst_idx);
nhist_ax(CC,'text','box','median','mode','noerror','color',[.8 .9 .8],'separate','axis',app.ax_burst_integral_his)
title(app.ax_burst_integral_his,'')

FF.burst_amplitude = app.burst_res.burst_amp(burst_idx);
nhist_ax(FF,'text','box','median','mode','noerror','color',[.8 .5 .8],'separate','axis',app.ax_burst_amplitude_his)
title(app.ax_burst_amplitude_his,'')

DD.burst_duration=(app.burst_res.burst_loc(burst_idx,2)-app.burst_res.burst_loc(burst_idx,1)) * app.burst_res.ts(2);
nhist_ax(DD,'text','box','median','mode','noerror','color',[.8 .8 1],'separate','axis',app.ax_burst_duration_his)
title(app.ax_burst_duration_his,'')

EE.inter_burst_interval = app.burst_res.dt_burst(burst_idx(1:end-1));
nhist_ax(EE,'text','box','median','mode','noerror','color',[.9 .8 .8],'separate','axis',app.ax_burst_interval_his)
title(app.ax_burst_interval_his,'')

app.lbl_num_burst.Text = ['Number of Bursts: ' num2str(sum(burst_idx)) ];
app.lbl_burst_rate.Text = ['Rate: ' num2str(sum(burst_idx)/dur) ' [Hz]'];
if ~isempty(app.hb_res)
    app.lbl_incidence.Text = ['Incidence: ' num2str(sum(burst_idx)/sum(hb_idx)) ' bursts/beats'];
else
    app.lbl_incidence.Text = 'no ECG available';
end
app.lbl_threshold.Text = ['HThreshold: ' num2str(app.burst_res.htresh) ' [' char(string(app.data(app.settings.channel_idx.msna).unit)) ']'];

end


