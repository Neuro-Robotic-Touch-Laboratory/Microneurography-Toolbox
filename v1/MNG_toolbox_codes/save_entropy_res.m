function save_entropy_res(app)
%save_entropy_res Summary of this function goes here
%   Detailed explanation goes here
[int_idxs,~] = listdlg('PromptString',{'Please select intervals ',...
    'to be plotted/saved.',''},...
    'SelectionMode','multiple','ListString',app.popup_int_entropy.Items);

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg','.eps'});
%path = uigetdir;
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';

[data,ts,name, ~] = current_signal(app, app.entropy_res(1).channel_idx);
data(:,2) = ts(1):ts(1):ts(2);
name = char(name);
tmp =find(vertcat(app.burst_ints.type) == 1);
int_names = {'full interval'};
borders = nan(length(tmp)+1,2);
borders(1,:) = [1 size(data,1)];

for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    
    borders(i+1,:) = [ceil(borders(i+1,1)/ts(1)), floor(borders(i+1,2)/ts(1))];
    int_names{i+1} = app.burst_ints(tmp(i)).name;
end

for j = 1:length(int_idxs)
    h = figure('Visible','off');
    set(h, 'NumberTitle', 'off', 'Name', int_names{int_idxs(j)});
    
    data_int = data(borders(int_idxs(j),1):borders(int_idxs(j),2),:);
    switch app.entropy_res(1).cfg.method
        case 'PE'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of permutation entropy ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int(:,2), data_int(:,1)],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata) + 1:end, 2 ), app.entropy_res(int_idxs(j)).outdata']};
        case 'PEeq'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of permutation entropy for ordinal patterns with tied ranks ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int(:,2), data_int(:,1)],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata) + 1:end, 2 ), app.entropy_res(int_idxs(j)).outdata']};
        case 'CE'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of conditional entropy of ordinal patterns ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int(:,2), data_int(:,1)]...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata) + 1:end, 2 ), app.entropy_res(int_idxs(j)).outdata']};
        case 'rePE'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of robust permutation entropy ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int(:,2), data_int(:,1)],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata) + 1:end, 2 ), app.entropy_res(int_idxs(j)).outdata']};
        case 'opdPE'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of ordinary pattern distribution permutation entropy ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int(:,2), data_int(:,1)],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata.ePe) + 1:end, 2 ), app.entropy_res(int_idxs(j)).outdata.ePe']};
        case 'all'
            title_str ={['Original time series ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of permutation entropy ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of conditional entropy of ordinal patterns ' name ' ' int_names{int_idxs(j)}],...
                        ['Values of robust permutation entropy ' name ' ' int_names{int_idxs(j)}]};
            plot_data = {[data_int( end - length( app.entropy_res(int_idxs(j)).outdata.PE ) + 1:end, 2), data_int( end - length( app.entropy_res(int_idxs(j)).outdata.PE ) + 1:end, 1)],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata.PE ) + 1:end, 2), app.entropy_res(int_idxs(j)).outdata.PE'],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata.CE ) + 1:end, 2), app.entropy_res(int_idxs(j)).outdata.CE'],...
                         [data_int( end - length( app.entropy_res(int_idxs(j)).outdata.RE ) + 1:end, 2), app.entropy_res(int_idxs(j)).outdata.RE']};
    
    end
    
    for i = 1: length( title_str)
        subplot(length(title_str),1,i)
        plot (plot_data{i}(:,1),plot_data{i}(:,2), 'k', 'LineWidth', .2 )
        title(title_str{i})
        xlabel('seconds')
    end
    
    for i = 1 : length(form_idxs)
        switch form_idxs(i)
            case 1
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_entropy_analysis_' simple_name(int_names{int_idxs(j)}) '_SIG_' name '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_entropy_analysis_' simple_name(int_names{int_idxs(j)}) '_SIG_' name '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_entropy_analysis_' simple_name(int_names{int_idxs(j)}) '_SIG_' name '.epsc'])
        end
    end
    close(h)
    app.lbl_working.Text = [num2str(round((100*j)/length(int_idxs))) '% done'];
end
end

