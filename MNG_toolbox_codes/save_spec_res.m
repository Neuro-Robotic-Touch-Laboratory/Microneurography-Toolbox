function save_spec_res(app)
%UNTITLED save spectral analysis results
%   Detailed explanation goes here

cols = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980];...
        [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560];...
        [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};

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

channel_idx = find(strcmp(app.popup_signal_spect.Value,app.popup_signal_spect.Items));
[data,ts,name, unit] = current_signal(app, channel_idx);
data(:,2) = ts(1):ts(1):ts(2);

int_names = {'full interval'};
borders = nan(length(tmp)+1,2);
borders(1,:) = [1 size(data,1)];
tmp =find(vertcat(app.burst_ints.type) == 1);
for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    
    borders(i+1,:) = [ceil(borders(i+1,1)/ts(1)), floor(borders(i+1,2)/ts(1))];
    int_names{i+1} = app.burst_ints(tmp(i)).name;
end


for i = 1:length(int_idxs)
    h = figure('Position', get(0, 'Screensize'),'Visible','on');
    set(h, 'NumberTitle', 'off', ...
    'Name', int_names{int_idxs(i)});
    
    subplot (5,2,1:2)
    hold on
    labels = {['DS = ' num2str(app.edt_ds1.Value)];['DS = ' num2str(app.edt_ds2.Value)]};
    title([name ' - ' int_names{1,i}])
    ylabel(unit)
    xlabel('time [s]')
%     for j = [app.edt_ds1.Value, app.edt_ds2.Value]
%         plot (downsample(data(borders(1):borders(2),2),j),downsample(data(borders(1):borders(2),1),j))
%     end
    ds = [app.edt_ds1.Value, app.edt_ds2.Value];
    for j = 1:2 
        plot (downsample(data(borders(int_idxs(i),1):borders(int_idxs(i),2),2),ds(j)),downsample(data(borders(int_idxs(i),1):borders(int_idxs(i),2),1),ds(j)),'Color',cols{1,j})
    end
    hold off 
    legend (labels)
    
    leg_str_1_1 ={[]};
    leg_str_2_1 ={[]};
    leg_str_3_1 ={[]};
    leg_str_4_1 ={[]};
    leg_str_1_2 ={[],[],[],[],[],[]};
    leg_str_2_2 ={[],[],[],[],[],[]};
    leg_str_3_2 ={[],[],[],[],[],[]};
    leg_str_4_2 ={[],[],[],[],[],[]};
    rem_idx = [];
    mx = nan(4,1);
    for j = 1 : 2
        subplot(5,2,3)
        hold on
        tmp_dt = mean(diff(app.spec_res.x_1_1{j,int_idxs(i)}));
        lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
        plot(app.spec_res.x_1_1{j,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_1_1{j,int_idxs(i)}(lim(1):lim(2)))
        mx(1) = max([mx(1),max(app.spec_res.y_1_1{j,int_idxs(i)}(int16(0.4/tmp_dt):int16(app.edt_max_frq.Value/tmp_dt)))]);
        leg_str_1_1{j} = char(string(app.spec_res.lbl_1_1{j,int_idxs(i)}));
        if j ==2
            legend (leg_str_1_1) 
        end
        hold off

        subplot(5,2,5)
        hold on
        tmp_dt = mean(diff(app.spec_res.x_2_1{j,int_idxs(i)}));
        lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
        plot(app.spec_res.x_2_1{j,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_2_1{j,int_idxs(i)}(lim(1):lim(2)))
        mx(2) = max([mx(2),max(app.spec_res.y_2_1{j,int_idxs(i)}(int16(0.3/tmp_dt):int16(app.edt_max_frq.Value/tmp_dt)))]);
        leg_str_2_1{j} = char(string(app.spec_res.lbl_2_1{j,int_idxs(i)}));
        if j ==2
            legend (leg_str_2_1) 
        end
        hold off

        subplot(5,2,7)
        hold on
        tmp_dt = mean(diff(app.spec_res.x_3_1{j,int_idxs(i)}));
        lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
        plot(app.spec_res.x_3_1{j,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_3_1{j,int_idxs(i)}(lim(1):lim(2)))
        mx(3) = max([mx(3),max(app.spec_res.y_3_1{j,int_idxs(i)}(int16(0.3/tmp_dt):int16(app.edt_max_frq.Value/tmp_dt)))]);
        leg_str_3_1{j} = char(string(app.spec_res.lbl_3_1{j,int_idxs(i)}));
        if j ==2
            legend (leg_str_3_1) 
        end
        hold off

        subplot(5,2,9)
        hold on
        tmp_dt = mean(diff(app.spec_res.x_4_1{j,int_idxs(i)}));
        lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
        plot(app.spec_res.x_4_1{j,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_4_1{j,int_idxs(i)}(lim(1):lim(2)))
        mx(4) = max([mx(4),max(app.spec_res.y_4_1{j,int_idxs(i)}(int16(0.3/tmp_dt):int16(app.edt_max_frq.Value/tmp_dt)))]);
        leg_str_4_1{j} = char(string(app.spec_res.lbl_4_1{j,int_idxs(i)}));
        if j ==2
            legend (leg_str_4_1) 
        end
        hold off
        
        for k = 1:3
            switch k
                case 1
                    do_plot = app.chkbx_show_o1.Value;
                case 2
                    do_plot = app.chkbx_show_o2.Value;
                case 3
                    do_plot = app.chkbx_show_o3.Value;
            end
            if do_plot
            
                plot_idx = k+(j-1)*3;
    
                subplot(5,2,4)
                hold on
                tmp_dt = mean(diff(app.spec_res.x_1_2{plot_idx,int_idxs(i)}));
                lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
                plot(app.spec_res.x_1_2{plot_idx,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_1_2{plot_idx,int_idxs(i)}(lim(1):lim(2)))
                leg_str_1_2{plot_idx} = char(string(app.spec_res.lbl_1_2{plot_idx,int_idxs(i)}));
                if j ==2
                    leg_str_1_2(rem_idx) = [];
                    legend (leg_str_1_2) 
                end
                hold off
                
                subplot(5,2,6)
                hold on
                tmp_dt = mean(diff(app.spec_res.x_2_2{plot_idx,int_idxs(i)}));
                lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
                plot(app.spec_res.x_2_2{plot_idx,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_2_2{plot_idx,int_idxs(i)}(lim(1):lim(2)))
                leg_str_2_2{plot_idx} = char(string(app.spec_res.lbl_2_2{plot_idx,int_idxs(i)}));
                if j ==2
                    leg_str_2_2(rem_idx) = [];
                    legend (leg_str_2_2) 
                end
                hold off
    
                subplot(5,2,8)
                hold on 
                tmp_dt = mean(diff(app.spec_res.x_3_2{plot_idx,int_idxs(i)}));
                lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
                plot(app.spec_res.x_3_2{plot_idx,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_3_2{plot_idx,int_idxs(i)}(lim(1):lim(2)))
                leg_str_3_2{plot_idx} = char(string(app.spec_res.lbl_3_2{plot_idx,int_idxs(i)}));
                if j ==2
                    leg_str_3_2(rem_idx) = [];
                    legend (leg_str_3_2) 
                end
                hold off
    
                subplot(5,2,10)
                hold on 
                tmp_dt = mean(diff(app.spec_res.x_4_2{plot_idx,int_idxs(i)}));
                lim= [round(app.edt_min_frq.Value/tmp_dt)  round(app.edt_max_frq.Value/tmp_dt)];
                plot(app.spec_res.x_4_2{plot_idx,int_idxs(i)}(lim(1):lim(2)), app.spec_res.y_4_2{plot_idx,int_idxs(i)}(lim(1):lim(2)))
                leg_str_4_2{plot_idx} = char(string(app.spec_res.lbl_4_2{plot_idx,int_idxs(i)}));
                if j ==2
                    leg_str_4_2(rem_idx) = [];
                    legend (leg_str_4_2) 
                end
                hold off
            else
                disp 'else case' 
                if j == 1
                    disp 'if executetd'
                    rem_idx = [rem_idx, k,k+3];
                end
            end
        end
    end
    lim= [app.edt_min_frq.Value  app.edt_max_frq.Value];
    subplot(5,2,3)
    %ym =1.5*max( app.spec_res.y_1_1{j,i}( int16(0.2/mean(diff(app.spec_res.x_1_1{j,i}))) : int16(app.edt_max_frq.Value/mean(diff(app.spec_res.x_1_1{j,i})))));
    ylim([0, 1.5*mx(1)])
    xlim(lim)
    subplot(5,2,5)
    %ym = 1.5*max( app.spec_res.y_2_1{j,i}( int16(0.2/mean(diff(app.spec_res.x_2_1{j,i}))) : int16(app.edt_max_frq.Value/mean(diff(app.spec_res.x_2_1{j,i})))));
    ylim([0, 1.5*mx(2)])
    xlim(lim)
    subplot(5,2,7)
    %ym = 1.5*max( app.spec_res.y_3_1{j,i}(int16(0.2/mean(diff(app.spec_res.x_3_1{j,i}))) : int16(app.edt_max_frq.Value/mean(diff(app.spec_res.x_3_1{j,i})))));
    ylim([0, 1.5*mx(3)])
    xlim(lim)
    subplot(5,2,9)
    %ym = 1.5*max( app.spec_res.y_4_1{j,i}( int16(0.2/mean(diff(app.spec_res.x_4_1{j,i}))) : int16(app.edt_max_frq.Value/mean(diff(app.spec_res.x_4_1{j,i})))));
    ylim([0, 1.5*mx(4)])
    xlim(lim)
    subplot(5,2,4)
    ylim([0 app.edt_amp.Value])
    xlim(lim)
    subplot(5,2,6)
    ylim([0 app.edt_amp.Value])
    xlim(lim)
    subplot(5,2,8)
    ylim([0 app.edt_amp.Value])
    xlim(lim)
    subplot(5,2,10)
    ylim([0 app.edt_amp.Value])
    xlim(lim)

    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_spectral_analysis_' simple_name(int_names{int_idxs(i)}) '_SIG_' name '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_spectral_analysis_' simple_name(int_names{int_idxs(i)}) '_SIG_' name '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_spectral_analysis_' simple_name(int_names{int_idxs(i)}) '_SIG_' name '.epsc'])
        end
    end
    app.lbl_working.Text = [ num2str(round(100*i /length(int_idxs))) '% done'];
end

end

