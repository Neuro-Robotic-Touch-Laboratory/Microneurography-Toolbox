function update_correlation_axis(app)

if isempty(app.corre_res)
    yyaxis(app.ax_signals_corre,'left')
    cla(app.ax_signals_corre)
    yyaxis(app.ax_signals_corre,'right')
    cla(app.ax_signals_corre)
    cla(app.ax_signal1_corre)
    cla(app.ax_signal2_corre)
    cla(app.ax_corre_1_1_1)
    cla(app.ax_corre_1_1_2)
    cla(app.ax_corre_1_2_1)
    cla(app.ax_corre_1_2_2)
    cla(app.ax_corre_2_1_1)
    cla(app.ax_corre_2_1_2)
    cla(app.ax_corre_2_2_1)
    cla(app.ax_corre_2_2_2)
    cla(app.ax_corre_3_1_1)
    cla(app.ax_corre_3_1_2)
    cla(app.ax_corre_3_2_1)
    cla(app.ax_corre_3_2_2)
else

    int_idx = find(strcmp(app.popup_int_corre.Value, app.popup_int_corre.Items));
    if int_idx ~=1
        tmp =find(vertcat(app.burst_ints.type) == 1);
        border = app.burst_ints(tmp(int_idx-1)).borders;
        int_name = app.burst_ints(tmp(int_idx-1)).name;    
    else
        border = [nan, nan];
        int_names = {'full interval'};
    end
    
    ds = app.edt_downsampling.Value;
    
    ch1_idx = find(strcmp(app.popup_coorelation_signal1.Value, app.popup_coorelation_signal1.Items));
    [data_1,ts_1,name_1, unit_1] = current_signal(app, ch1_idx);
    data_1(:,2) = ts_1(1):ts_1(1):ts_1(2);
    
    ch2_idx = find(strcmp(app.popup_coorelation_signal2.Value, app.popup_coorelation_signal2.Items));
    [data_2,ts_2,name_2, unit_2] = current_signal(app, ch2_idx);
    data_2(:,2) = ts_2(1):ts_2(1):ts_2(2);
    
    if ts_1(1) > ts_2(1)
        idx = find(data_2(:,2) == data_1(1,2));
        data_2(1:idx-1,:) = [];
    else
        idx = find(data_1(:,2) == data_2(1,2));
        data_1(1:idx-1,:) = [];
    end
    
    data_1_ds = [downsample(data_1(:,1), .01/ts_1(1)),downsample(data_1(:,2), .01/ts_1(1))];
    
    data_2_ds = [downsample(data_2(:,1), .01/ts_2(1)),downsample(data_2(:,2), .01/ts_2(1))];
    
    % data_2_ds(:,2) = data_1_ds(:,2);
    
    data_1_dsds = [downsample(data_1_ds(:,1), ds), downsample(data_1_ds(:,2), ds)];
    data_2_dsds = [downsample(data_2_ds(:,1), ds), downsample(data_2_ds(:,2), ds)];

    if app.chkbx_movmean.Value
        data_1_dsds(:,1) = movmean(data_1_dsds(:,1),app.edt_movmean.Value /(.01*ds));
        data_2_dsds(:,1) = movmean(data_2_dsds(:,1),app.edt_movmean.Value /(.01*ds));
    end

    border(border(:,1) < 1 ,1) = 1;
    border(border(:,2) > floor(min([data_1_dsds(end,2), data_2_dsds(end,2)])),2) = floor(min([data_1_dsds(end,2), data_2_dsds(end,2)]));
    
    if int_idx == 1
        data_1_dsds_int = data_1_dsds;
        data_2_dsds_int = data_2_dsds;
    else
        data_1_dsds_int = data_1_dsds(int32(border(1)*1/(.01*ds)):int32(border(2)*1/(.01*ds)),:);
        data_2_dsds_int = data_2_dsds(int32(border(1)*1/(.01*ds)):int32(border(2)*1/(.01*ds)),:);
    end

    yyaxis(app.ax_signals_corre,'left')
    plot(app.ax_signals_corre, data_1_dsds_int(:,2), data_1_dsds_int(:,1),'Linewidth',2 )
    %ylabel([name_1 ' ' unit_1]),xlabel('time [s]')
    %hold(app.ax_signals_corre, 'on')
    yyaxis(app.ax_signals_corre, 'right') 
    plot(app.ax_signals_corre, data_2_dsds_int(:,2), data_2_dsds_int(:,1),'Linewidth',2 )
    %ylabel([name_2 ' ' unit_2]))
    tlim=[min(data_1_dsds_int(1,2),data_2_dsds_int(1,2)) max(data_1_dsds_int(end,2),data_2_dsds_int(end,2))];
    xlim(app.ax_signals_corre, tlim);
    %legend(name_1, name_2)
    %title([name_1 ' & ' name2])
    app.lbl_signals_corre1.Text = name_1;
    app.lbl_signals_corre2.Text = name_2;
    app.lbl_signals_corre.Text = ' & ';


    plot(app.ax_signal1_corre, data_1_dsds_int(:,2), data_1_dsds_int(:,1),'.','MarkerSize',12);
    hold(app.ax_signal1_corre, 'on')
    plot(app.ax_signal1_corre, data_1_dsds_int(:,2), data_1_dsds_int(:,1),'-','LineWidth',0.5,'color',[0.3 0.3 0.3]);
    %title(name_1);
    app.lbl_signal1_corre.Text = name_1; 
    xlim(app.ax_signal1_corre, tlim)
    hold(app.ax_signal1_corre, 'off')
%     colorbar('off')
    %xlabel('time [s]')
    % ylabel('Period [s]')


    plot(app.ax_signal2_corre, data_2_dsds_int(:,2), data_2_dsds_int(:,1),'.','color',[0.8500 0.3250 0.0980],'MarkerSize',12);
    hold(app.ax_signal2_corre, 'on')
    plot(app.ax_signal2_corre, data_2_dsds_int(:,2), data_2_dsds_int(:,1),'-','LineWidth',0.5,'color',[0.3 0.3 0.3]);
    %title(name_2)
    app.lbl_signal2_corre.Text = name_2;
    xlim(app.ax_signal2_corre, tlim)
    hold(app.ax_signal2_corre, 'off')
%     colorbar('off')
%     xlabel('time [s]')
    % ylabel('Period [s]')

    lagdata=[app.corre_res(int_idx).lag1, app.corre_res(int_idx).lag2, app.corre_res(int_idx).lag3];
    axs(:,:,1) = [app.ax_corre_1_1_1, app.ax_corre_1_1_2;...
                  app.ax_corre_1_2_1, app.ax_corre_1_2_2]; 
    axs(:,:,2) = [app.ax_corre_2_1_1, app.ax_corre_2_1_2;...
                  app.ax_corre_2_2_1, app.ax_corre_2_2_2];
    axs(:,:,3) = [app.ax_corre_3_1_1, app.ax_corre_3_1_2;...
                  app.ax_corre_3_2_1, app.ax_corre_3_2_2];
    lbls = [app.lbl_corre_1_1_2, app.lbl_corre_1_2_1, app.lbl_xcorre1;...
            app.lbl_corre_2_1_2, app.lbl_corre_2_2_1, app.lbl_xcorre2;...
            app.lbl_corre_3_1_2, app.lbl_corre_3_2_1, app.lbl_xcorre3];
    lags = [app.edt_corre_lag1.Value, app.edt_corre_lag2.Value, app.edt_corre_lag3.Value];

    for i = 1:3
        lbls(i,3).Text = ['Lag' num2str(i) ': ' num2str(lags(i))];
        histogram(axs(1,1,i), lagdata(i).ax11.Data ,lagdata(i).ax11.BinEdges,'FaceColor',[0,0,1],'EdgeColor',[0,1,1])
        xlim(axs(1,1,i), lagdata(i).ax11.XLim)
        ylim(axs(1,1,i), lagdata(i).ax11.YLim)
        set(axs(1,1,i),'xtick',[],'ytick',[],'xgrid','off','ygrid','off','Toolbar',[]);
        
        plot(axs(1,2,i), lagdata(i).ax12.line2.XData, lagdata(i).ax12.line2.YData, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'none')
        set(axs(1,2,i),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off','xtick',[]);
        hold(axs(1,2,i), 'on')
        plot(axs(1,2,i), lagdata(i).ax12.line1.XData, lagdata(i).ax12.line1.YData,'-m', 'LineWidth', 0.5000)
        lbls(i,1).Text = lagdata(i).ax12.text.String;
        xlim(axs(1,2,i), lagdata(i).ax12.XLim)
        ylim(axs(1,2,i), lagdata(i).ax12.YLim)
        set(axs(1,2,i),'xtick',[],'ytick',[],'xgrid','off','ygrid','off','Toolbar',[]);
        hold(axs(1,2,i), 'off')

        plot(axs(2,1,i), lagdata(i).ax21.line2.XData, lagdata(i).ax21.line2.YData, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'none')
        set(axs(2,1,i),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off','xtick',[]);
        hold(axs(2,1,i), 'on')
        plot(axs(2,1,i), lagdata(i).ax21.line1.XData, lagdata(i).ax21.line1.YData,'-m', 'LineWidth', 0.5000)
        lbls(i,2).Text = lagdata(i).ax21.text.String;
        xlim(axs(2,1,i), lagdata(i).ax21.XLim)
        ylim(axs(2,1,i), lagdata(i).ax21.YLim)
        set(axs(2,1,i),'xtick',[],'ytick',[],'xgrid','off','ygrid','off','Toolbar',[]);
        hold(axs(2,1,i), 'off')
        
        histogram(axs(2,2,i), lagdata(i).ax22.Data ,lagdata(i).ax22.BinEdges,'FaceColor',[0,0,1],'EdgeColor',[0,1,1])
        xlim(axs(2,2,i), lagdata(i).ax22.XLim)
        ylim(axs(2,2,i), lagdata(i).ax22.YLim)
        set(axs(2,2,i),'xtick',[],'ytick',[],'xgrid','off','ygrid','off','Toolbar',[]);

    end

   
end
end