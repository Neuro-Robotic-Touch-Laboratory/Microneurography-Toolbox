int_idx = 1;
xl = [app.burst_res.ts(1), app.burst_res.ts(3)];
figure 
ax(1) = subplot (2,1,1);



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
hold off 


plot(t_msna(2001:end), app.burst_res.x(2001:end),'LineWidth',1.2,'HitTest','off')  
ylim(yl)  
hold('on')

hb = nan(length(app.hb_res.t_events(hb_idx))*3,2);
tmp_idx = 1:3:length(app.hb_res.t_events(hb_idx))*3;
hb(tmp_idx,1)= app.hb_res.t_events(hb_idx);
hb(tmp_idx+1,1)= app.hb_res.t_events(hb_idx);
yl = ylim(app.ax_msna_int);
hb(tmp_idx,2)= yl(1);
hb(tmp_idx+1,2)= yl(2);
plot (hb(:,1),hb(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2],'HitTest','off')

xlim(xl)

hold('on')
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

plot(burst_plot1(:,1),burst_plot1(:,2),'Color',colors(1,:),'LineWidth',1.5)
plot(burst_plot2(:,1),burst_plot2(:,2),'Color',colors(2,:),'LineWidth',1.5)
plot(burst_plot3(:,1),burst_plot3(:,2),'Color',colors(3,:),'LineWidth',1.5)
plot(burst_plot4(:,1),burst_plot4(:,2),'Color',colors(4,:),'LineWidth',1.5)
plot(burst_plot5(:,1),burst_plot5(:,2),'Color',colors(5,:),'LineWidth',1.5)
plot(burst_plot_rem(:,1),burst_plot_rem(:,2),'Color',[.7,.7,.7],'LineWidth',1.5)


ax(2) = subplot (2,1,2);


hb_idx = app.hb_res.use_beats(:,1) ;

plot (t_msna,x)
hold on 

hb = nan(length(app.hb_res.t_events(hb_idx))*3,2);
tmp_idx = 1:3:length(app.hb_res.t_events(hb_idx))*3;
hb(tmp_idx,1)= app.hb_res.t_events(hb_idx);
hb(tmp_idx+1,1)= app.hb_res.t_events(hb_idx);
yl = ylim(app.ax_msna_int);
hb(tmp_idx,2)= yl(1);
hb(tmp_idx+1,2)= yl(2);
plot (hb(:,1),hb(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2],'HitTest','off')

        
colors = [1.0000, 0, 0;...
          0, 1.0000, 0;...
          0, 0, 0.1724;...
          1.0000, 0.1034, 0.7241;...
          1.0000, 0.8276, 0];

burst_plot1 =[nan,nan]; burst_plot2 =[nan,nan]; burst_plot3 =[nan,nan]; burst_plot4 =[nan,nan]; burst_plot5 =[nan,nan]; burst_plot_rem = [nan,nan];
burst_amp_plot1 = [nan,nan]; burst_amp_plot2 = [nan,nan]; burst_amp_plot3 = [nan,nan]; burst_amp_plot4 = [nan,nan]; burst_amp_plot5 = [nan,nan]; burst_amp_plot_rem = [nan,nan]; 
tmp = find(brsts(:,8) == 1);
tmp_brsts = brsts(tmp,:);

for j = 1:5: ceil(size(tmp_brsts,1)/5)*5 %%%% floor
    if j <= size(tmp_brsts,1)
        burst_plot1 = [burst_plot1;[(t_msna(tmp_brsts(j,1):tmp_brsts(j,3)))',...
                       (x(tmp_brsts(j,1):tmp_brsts(j,3)))'];...
                       [t_msna(tmp_brsts(j,1)), x(tmp_brsts(j,1))];...
                       [nan,nan]];
    end
    if j+1 <= size(tmp_brsts,1)  
        burst_plot2 = [burst_plot2;[(t_msna(tmp_brsts(j+1,1):tmp_brsts(j+1,3)))',...
                       (x(tmp_brsts(j+1,1):tmp_brsts(j+1,3)))'];...
                       [t_msna(tmp_brsts(j+1,1)), x(tmp_brsts(j+1,1))];...
                       [nan,nan]];
    end
    if j+2 <= size(tmp_brsts,1)
        burst_plot3 = [burst_plot3;[(t_msna(tmp_brsts(j+2,1):tmp_brsts(j+2,3)))',...
                       (x(tmp_brsts(j+2,1):tmp_brsts(j+2,3)))'];...
                       [t_msna(tmp_brsts(j+2,1)), x(tmp_brsts(j+2,1))];...
                       [nan,nan]];
    end
    if j+3 <= size(tmp_brsts,1)
        burst_plot4 = [burst_plot4;[(t_msna(tmp_brsts(j+3,1):tmp_brsts(j+3,3)))',...
                       (x(tmp_brsts(j+3,1):tmp_brsts(j+3,3)))'];...
                       [t_msna(tmp_brsts(j+3,1)), x(tmp_brsts(j+3,1))];...
                       [nan,nan]];
    end
    if j+4 <= size(tmp_brsts,1) 
        burst_plot5 = [burst_plot5;[(t_msna(tmp_brsts(j+4,1):tmp_brsts(j+4,3)))',...
                       (x(tmp_brsts(j+4,1):tmp_brsts(j+4,3)))'];...
                       [t_msna(tmp_brsts(j+4,1)), x(tmp_brsts(j+4,1))];...
                       [nan,nan]];
    end
end
tmp = find(brsts(:,8)==0);
tmp_brsts = brsts(tmp,:);

for j = 1: length(tmp)
    burst_plot_rem = [burst_plot_rem;[(t_msna(tmp_brsts(j,1):tmp_brsts(j,3)))',...
                       (x(tmp_brsts(j,1):tmp_brsts(j,3)))'];...
                       [t_msna(tmp_brsts(j,1)), x(tmp_brsts(j,1))];...
                       [nan,nan]];
end
hold on 
plot(burst_plot1(:,1),burst_plot1(:,2),'Color',colors(1,:),'LineWidth',1.5)
plot(burst_plot2(:,1),burst_plot2(:,2),'Color',colors(2,:),'LineWidth',1.5)
plot(burst_plot3(:,1),burst_plot3(:,2),'Color',colors(3,:),'LineWidth',1.5)
plot(burst_plot4(:,1),burst_plot4(:,2),'Color',colors(4,:),'LineWidth',1.5)
plot(burst_plot5(:,1),burst_plot5(:,2),'Color',colors(5,:),'LineWidth',1.5)
plot(burst_plot_rem(:,1),burst_plot_rem(:,2),'Color',[.7,.7,.7],'LineWidth',1.5)


linkaxes(ax,'x')
