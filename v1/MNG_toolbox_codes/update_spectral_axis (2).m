function update_spectral_axis(app,static)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
int_idx = find(strcmp(app.popup_int_spect.Value,app.popup_int_spect.Items));
% cols = {[0 0 1],[1 0 0];[0 1 0],[0 1 1]; [1 0 1],[1 1 0]};
cols = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980];...
        [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560];...
        [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};

if static
    
    labels = {['DS = ' num2str(app.edt_ds1.Value)];['DS = ' num2str(app.edt_ds2.Value)]};
    
    channel_idx = find(strcmp(app.popup_signal_spect.Value,app.popup_signal_spect.Items));
    
   
    cla(app.ax_signal)
    [data,ts,name, unit] = current_signal(app, channel_idx);
    data(:,2) = ts(1):ts(1):ts(2);
    if int_idx == 1
        borders = [1 size(data,1)];
    else
        
        tmp =find(vertcat(app.burst_ints.type) == 1);
        borders = app.burst_ints(tmp(int_idx-1)).borders;
        borders = [ceil(borders(1)/ts(1)), floor(borders(2)/ts(1))];
    end
    title(app.ax_signal, name)
    ylabel(app.ax_signal,unit)
    hold(app.ax_signal,'on')
    ds = [app.edt_ds1.Value, app.edt_ds2.Value];
    for i = 1:2 
        plot (app.ax_signal,downsample(data(borders(1):borders(2),2),ds(i)),downsample(data(borders(1):borders(2),1),ds(i)),'Color',cols{1,i})
    end
    hold(app.ax_signal, 'off')
    legend (app.ax_signal, labels)
end
cla(app.ax_spect_1_1)
cla(app.ax_spect_2_1)
cla(app.ax_spect_3_1)
cla(app.ax_spect_4_1)
cla(app.ax_spect_1_2)
cla(app.ax_spect_2_2)
cla(app.ax_spect_3_2)
cla(app.ax_spect_4_2)
hold(app.ax_spect_1_1, 'on')
hold(app.ax_spect_2_1, 'on')
hold(app.ax_spect_3_1, 'on')
hold(app.ax_spect_4_1, 'on')
hold(app.ax_spect_1_2, 'on')
hold(app.ax_spect_2_2, 'on')
hold(app.ax_spect_3_2, 'on')
hold(app.ax_spect_4_2, 'on')

if ~isempty(app.spec_res)

    for i = 1:2
        xl = [app.edt_min_frq.Value, app.edt_max_frq.Value];
        xlim(app.ax_spect_1_1, xl)
%         tmp_dt = mean(diff(app.spec_res.x_1_1{i,int_idx}));
%         lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
%         if lim(1)<= 1
%             lim(1) = 1;
%         end
        plot(app.ax_spect_1_1,app.spec_res.x_1_1{i,int_idx},...     %(lim(1):lim(2)), ...
                              app.spec_res.y_1_1{i,int_idx},'Color',cols{1,i})        %(lim(1):lim(2)))
        
        ym =1.5*max( app.spec_res.y_1_1{i,int_idx}( int16(0.2/mean(diff(app.spec_res.x_1_1{i,int_idx}))):end));
        ylim(app.ax_spect_1_1, [0, ym])

%         tmp_dt = mean(diff(app.spec_res.x_2_1{i,int_idx}));
%         lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
%         if lim(1)<= 1
%             lim(1) = 1;
%         end
        plot(app.ax_spect_2_1,app.spec_res.x_2_1{i,int_idx},...     %(lim(1):lim(2)),...
                              app.spec_res.y_2_1{i,int_idx},'Color',cols{1,i})        %(lim(1):lim(2)))
%         xlim(app.ax_spect_2_1, xl)
        ym = 1.5*max( app.spec_res.y_2_1{i,int_idx}( int16(0.2/mean(diff(app.spec_res.x_2_1{i,int_idx}))):end));
        ylim(app.ax_spect_2_1, [0, ym])

%         tmp_dt = mean(diff(app.spec_res.x_3_1{i,int_idx}));
%         lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
%         if lim(1)<= 1
%             lim(1) = 1;
%         end
        plot(app.ax_spect_3_1,app.spec_res.x_3_1{i,int_idx},...     %(lim(1):lim(2)),...
                              app.spec_res.y_3_1{i,int_idx},'Color',cols{1,i})        %(lim(1):lim(2)))
%         xlim(app.ax_spect_3_1, xl)
        
        ym = 1.5*max( app.spec_res.y_3_1{i,int_idx}(0.2/mean(diff(app.spec_res.x_3_1{i,int_idx})):end));
        ylim(app.ax_spect_3_1, [0, ym])
%         tmp_dt = mean(diff(app.spec_res.x_4_1{i,int_idx}));
%         lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
%         if lim(1)<= 1
%             lim(1) = 1;
%         end
        plot(app.ax_spect_4_1,app.spec_res.x_4_1{i,int_idx},...     %(lim(1):lim(2)),...
                              app.spec_res.y_4_1{i,int_idx},'Color',cols{1,i})        %(lim(1):lim(2)))
        %xlim(app.ax_spect_4_1, xl)
        ym = 1.5*max( app.spec_res.y_4_1{i,int_idx}( int16(0.2/mean(diff(app.spec_res.x_4_1{i,int_idx}))):end));
        ylim(app.ax_spect_4_1, ...
                [0, ym])

        for j = 1:3
            switch j 
                case 1
                    do_plot = app.chkbx_show_o1.Value;
                case 2
                    do_plot = app.chkbx_show_o2.Value;
                case 3
                    do_plot = app.chkbx_show_o3.Value;
            end
            if do_plot
                plot_idx = j+(i-1)*3;
    %             tmp_dt = mean(diff(app.spec_res.x_1_2{plot_idx,int_idx}));
    %             lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
    %             if lim(1)<= 1
    %                 lim(1) = 1;
    %             end
                plot(app.ax_spect_1_2,app.spec_res.x_1_2{plot_idx,int_idx},...  %(lim(1):lim(2)),...
                                      app.spec_res.y_1_2{plot_idx,int_idx},'Color',cols{j,i})     %(lim(1):lim(2)))
                %xlim(app.ax_spect_1_2, xl)
                ylim(app.ax_spect_1_2, [0 app.edt_amp.Value])
    
    %             tmp_dt = mean(diff(app.spec_res.x_2_2{plot_idx,int_idx}));
    %             lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
    %             if lim(1)<= 1
    %                 lim(1) = 1;
    %             end
                plot(app.ax_spect_2_2,app.spec_res.x_2_2{plot_idx,int_idx},...  %(lim(1):lim(2)),...
                                      app.spec_res.y_2_2{plot_idx,int_idx},'Color',cols{j,i})     %(lim(1):lim(2)))
                %xlim(app.ax_spect_2_2, xl)
                ylim(app.ax_spect_2_2, [0 app.edt_amp.Value])
    
    %             tmp_dt = mean(diff(app.spec_res.x_3_2{plot_idx,int_idx}));
    %             lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
    %             if lim(1)<= 1
    %                 lim(1) = 1;
    %             end
                plot(app.ax_spect_3_2,app.spec_res.x_3_2{plot_idx,int_idx},...  %(lim(1):lim(2)), ...
                                  app.spec_res.y_3_2{plot_idx,int_idx},'Color',cols{j,i})         %(lim(1):lim(2)))
                %xlim(app.ax_spect_3_2, xl)
                ylim(app.ax_spect_3_2, [0 app.edt_amp.Value])
    
    %             tmp_dt = mean(diff(app.spec_res.x_4_2{plot_idx,int_idx}));
    %             lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
    %             if lim(1)<= 1
    %                 lim(1) = 1;
    %             end
                plot(app.ax_spect_4_2,app.spec_res.x_4_2{plot_idx,int_idx},...  %(lim(1):lim(2)), ...
                                      app.spec_res.y_4_2{plot_idx,int_idx},'Color',cols{j,i})     %(lim(1):lim(2)))
                %xlim(app.ax_spect_4_2, xl)
                ylim(app.ax_spect_4_2, [0 app.edt_amp.Value])
            end
        end
%         leg_str = {['ds: ' num2str(app.edt_ds1.Value)],...
%                    ['ds: ' num2str(app.edt_ds2.Value)]};
%         legend(app.ax_spect_1_1,leg_str)
%         legend(app.ax_spect_2_1,leg_str)
%         legend(app.ax_spect_3_1,leg_str)
%         legend(app.ax_spect_4_1,leg_str)
        app.lbl_ds1_o1.Text = ['ds:' num2str(app.edt_ds1.Value) '; O:' num2str(app.edt_o1.Value)];
        app.lbl_ds1_o2.Text = ['ds:' num2str(app.edt_ds1.Value) '; O:' num2str(app.edt_o2.Value)];
        app.lbl_ds1_o3.Text = ['ds:' num2str(app.edt_ds1.Value) '; O:' num2str(app.edt_o3.Value)];
        app.lbl_ds2_o1.Text = ['ds:' num2str(app.edt_ds2.Value) '; O:' num2str(app.edt_o1.Value)];
        app.lbl_ds2_o2.Text = ['ds:' num2str(app.edt_ds2.Value) '; O:' num2str(app.edt_o2.Value)];
        app.lbl_ds2_o3.Text = ['ds:' num2str(app.edt_ds2.Value) '; O:' num2str(app.edt_o3.Value)];
        app.lbl_ds1.Text = ['ds:' num2str(app.edt_ds1.Value)];
        app.lbl_ds2.Text = ['ds:' num2str(app.edt_ds2.Value)];

    end   
end
hold(app.ax_spect_1_1, 'off')
hold(app.ax_spect_2_1, 'off')
hold(app.ax_spect_3_1, 'off')
hold(app.ax_spect_4_1, 'off')
hold(app.ax_spect_1_2, 'off')
hold(app.ax_spect_2_2, 'off')
hold(app.ax_spect_3_2, 'off')
hold(app.ax_spect_4_2, 'off')

end

