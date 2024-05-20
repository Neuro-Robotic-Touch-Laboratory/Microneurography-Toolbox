function corre_report(app)


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


if ~isempty(app.burst_ints)
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = nan(length(tmp)+1,2);
    int_names = {'full interval'};

    for i = 1: length(tmp)        
        borders(i+1,:) = app.burst_ints(tmp(i)).borders;
        
        int_names{i+1} = app.burst_ints(tmp(i)).name;
    end
else
    borders = [nan, nan];
    int_names = {'full interval'};
end

lags = [app.edt_corre_lag1.Value, app.edt_corre_lag2.Value, app.edt_corre_lag3.Value];
ds = app.edt_downsampling.Value;

ch1_idx = find(strcmp(app.popup_coorelation_signal1.Value, app.popup_coorelation_signal1.Items));
[data_1,ts_1,name_1, unit_1] = current_signal(app, ch1_idx);
data_1(:,2) = ts_1(1):ts_1(1):ts_1(2);
name_1 = char(string(name_1));

ch2_idx = find(strcmp(app.popup_coorelation_signal2.Value, app.popup_coorelation_signal2.Items));
[data_2,ts_2,name_2, unit_2] = current_signal(app, ch2_idx);
data_2(:,2) = ts_2(1):ts_2(1):ts_2(2);
name_2 = char(string(name_2));

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


for i = 1:length(int_idxs)
    
    if int_idxs(i) == 1
        data_1_dsds_int = data_1_dsds;
        data_2_dsds_int = data_2_dsds;
    else
        data_1_dsds_int = data_1_dsds(int32(borders(int_idxs(i),1)*1/(.01*ds)):int32(borders(int_idxs(i),2)*1/(.01*ds)),:);
        data_2_dsds_int = data_2_dsds(int32(borders(int_idxs(i),1)*1/(.01*ds)):int32(borders(int_idxs(i),2)*1/(.01*ds)),:);
    end

    h = figure('Position', get(0, 'Screensize'),'Visible','off');
    set(h, 'NumberTitle', 'off', ...
    'Name', int_names{int_idxs(i)});

    subplot(3,2,1)
    
    yyaxis left
    plot(data_1_dsds_int(:,2), data_1_dsds_int(:,1),'Linewidth',2 )
    ylabel([name_1 ' ' unit_1]),xlabel('time [s]')
    %hold(app.ax_signals_corre, 'on')
    yyaxis right
    plot(data_2_dsds_int(:,2), data_2_dsds_int(:,1),'Linewidth',2 )
    ylabel([name_2 ' ' unit_2])
    tlim=[min(data_1_dsds_int(1,2),data_2_dsds_int(1,2)) max(data_1_dsds_int(end,2),data_2_dsds_int(end,2))];
    xlim(tlim);
    legend(name_1, name_2)
    title([name_1 ' & ' name_2])
    
    subplot(3,2,3)

    plot(data_1_dsds_int(:,2), data_1_dsds_int(:,1),'.','MarkerSize',12);
    hold on
    plot(data_1_dsds_int(:,2), data_1_dsds_int(:,1),'-','LineWidth',0.5,'color',[0.3 0.3 0.3]);
    title(name_1); 
    xlim(tlim)
    hold off
    colorbar('off')
    xlabel('time [s]')
    ylabel(unit_1)
    
     subplot(3,2,5)

    plot(data_2_dsds_int(:,2), data_2_dsds_int(:,1),'.','color',[0.8500 0.3250 0.0980],'MarkerSize',12);
    hold on
    plot(data_2_dsds_int(:,2), data_2_dsds_int(:,1),'-','LineWidth',0.5,'color',[0.3 0.3 0.3]);
    title(name_2)
    xlim(tlim)
    hold off
    colorbar('off')
    xlabel('time [s]')
    ylabel(unit_2)

    subplot(3,2,2)
    lag01=lags(1)*1/(.01*ds);
    [~,~] = corrplot1([data_1_dsds_int(1:end-lag01,1),data_2_dsds_int(lag01+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    title(['Lag 1: ' num2str(app.edt_corre_lag1.Value) 's']) 

    
    subplot(3,2,4)
    lag02=lags(2)*1/(.01*ds);
    [~,~] = corrplot1([data_1_dsds_int(1:end-lag02,1),data_2_dsds_int(lag02+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    title(['Lag 2: ' num2str(app.edt_corre_lag2.Value) 's'])

    subplot(3,2,6)
    lag03=lags(3)*1/(.01*ds);
    [~,~] = corrplot1([data_1_dsds_int(1:end-lag03,1),data_2_dsds_int(lag03+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    title(['Lag 3: ' num2str(app.edt_corre_lag3.Value) 's'])

    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                %h.Visible = 'on';
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_correlation_' simple_name(int_names{int_idxs(i)}) '_SIG_' name_1 '+' name_2 '.fig'],'compact')
                %switch_vis([path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_correlation_' simple_name(int_names{int_idxs(i)}) '_SIG_' name_1 '+' name_2 '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_correlation_' simple_name(int_names{int_idxs(i)}) '_SIG_' name_1 '+' name_2 '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_correlation_' simple_name(int_names{int_idxs(i)}) '_SIG_' name_1 '+' name_2 '.epsc'])
        end
    end
    close(h)
    app.lbl_working.Text = [num2str(round((100*i)/length(int_idxs))) '% done'];
end