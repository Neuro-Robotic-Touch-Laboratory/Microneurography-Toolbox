function save_regression_res(app)

[int_idxs,~] = listdlg('PromptString',{'Please select intervals ',...
    'to be plotted/saved.',''},...
    'SelectionMode','multiple','ListString',app.popup_int_stepwise.Items);

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg','.eps'});
%path = uigetdir;
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';

for i = 1:length(int_idxs)
    if iscell(app.stepwise_res(int_idxs(i)).int_name)
        int_name = app.stepwise_res(int_idxs(i)).int_name{1,1};
    else
        int_name = app.stepwise_res(int_idxs(i)).int_name;
    end

    plot_cell = {[]};
    plot_cell{1,1} = 'par';
    plot_cell{1,2} = 'y estimated';
    tmp = num2cell( app.stepwise_res(int_idxs(i)).data.par );
    plot_cell(2:length(tmp)+1,1) = tmp; 
    tmp = num2cell( app.stepwise_res(int_idxs(i)).data.y_estimates );
    plot_cell(2:length(tmp)+1,2) = tmp; 
    writecell(plot_cell,...
              [path '\' file '_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_regression_results.xls' ],... 
              'FileType', 'spreadsheet', 'Sheet', int_name, 'Range', 'A1');
    
    plot_cell = {[]};
    plot_cell{1,1} = 'regressors';
    if iscell(app.stepwise_res(int_idxs(i)).data.reg)
        for j = 1:length(app.stepwise_res(int_idxs(i)).data.reg)
            plot_cell{1+j} = app.stepwise_res(int_idxs(i)).data.reg{1,j};
        end
    else
        plot_cell{2} = app.stepwise_res(int_idxs(i)).data.reg;
    end
    writecell(plot_cell,...
        [path '\' file '_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_regression_results.xls' ],...
        'FileType', 'spreadsheet', 'Sheet', int_name, 'Range', 'D1');

    plot_cell = {[]};
    plot_cell{1,1} = 'selected regressors';
    if iscell(app.stepwise_res(int_idxs(i)).data.sel_reg)
        for j = 1:length(app.stepwise_res(int_idxs(i)).data.sel_reg)
            plot_cell{1+j} = app.stepwise_res(int_idxs(i)).data.sel_reg{1,j};
        end
    else
        plot_cell{2} = app.stepwise_res(int_idxs(i)).data.sel_reg;
    end
    writecell(plot_cell,...
        [path '\' file '_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_regression_results.xls' ],...
        'FileType', 'spreadsheet', 'Sheet', int_name, 'Range', 'D2');

    plot_cell = {[]};
    plot_cell{1,1} = 'Beta' ;
    if iscell(app.stepwise_res(int_idxs(i)).data.beta)
        for j = 1:length(app.stepwise_res(int_idxs(i)).data.beta)
            plot_cell{1+j} = app.stepwise_res(int_idxs(i)).data.beta{1,j};
        end
    else
        plot_cell{2} = app.stepwise_res(int_idxs(i)).data.beta;
    end
    writecell(plot_cell,...
        [path '\' file '_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_regression_results.xls' ],...
        'FileType', 'spreadsheet', 'Sheet', int_name, 'Range', 'D3');
    
    h = figure('Position', get(0, 'Screensize'),'Visible', 'off');
    BlandAltman(h, app.stepwise_res(int_idxs(i)).data.par,  app.stepwise_res(int_idxs(i)).data.y_estimates)
  
    sgtitle([' Y-predicted: ' '____________' app.stepwise_res(int_idxs(i)).pred_name ' ' ' ' ' ' ' Regressors: ' '____________' app.stepwise_res(int_idxs(i)).data.sel_reg])
    
    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_stepwise_regression_' simple_name(app.stepwise_res(int_idxs(i)).int_name) '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_stepwise_regression_' simple_name(app.stepwise_res(int_idxs(i)).int_name) '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_stepwise_regression_' simple_name(app.stepwise_res(int_idxs(i)).int_name) '.epsc'])
        end
    end
    close(h)
    app.lbl_working.Text = [num2str(round((100*i)/length(int_idxs))) '% done'];
end 