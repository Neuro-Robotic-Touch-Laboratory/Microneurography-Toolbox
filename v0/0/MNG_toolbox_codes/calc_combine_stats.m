function calc_combine_stats(app)

intervals = {[]};
signals = {[]};


for i = 1: length(app.comb_stat)
    if i == 1
        intervals = app.comb_stat(i).intervals;
        signals = app.comb_stat(i).signals;
    else
        intervals = unique(horzcat(intervals, app.comb_stat(i).intervals),'stable');
        signals = unique(horzcat(signals, app.comb_stat(i).signals),'stable');
    end
end

collect_cell = {[]};
collect_cell{length(signals),length(intervals)} = [];
recording_cell = {[]};
recording_cell{length(signals),length(intervals)} = [];
recording_list = {[]};
recording_list{length(app.settings.brst_files)} = [];

for i = 1: length(app.settings.brst_files)
    tmp_cell = readcell(app.settings.brst_files(i).filename,'Sheet', 'general');
    recording_list{i} = [app.settings.brst_files(i).recording ' - ' app.settings.brst_files(i).interval];
    for j = 1 : length(signals)
        tmp_idx = find(contains(app.settings.brst_files(i).signals(:,1),signals{j}));
        if ~isempty(tmp_idx)
            sig_idx = app.settings.brst_files(i).signals{tmp_idx,2};
            
            for k = 1 : length(intervals)
                tmp_idx = find(contains(app.settings.brst_files(i).subintervals(:,1),intervals{k}));
                if ~isempty(tmp_idx)
                    int_idx = app.settings.brst_files(i).subintervals{tmp_idx,2};
                    if isempty(collect_cell{j,k})
                        collect_cell{j,k} = tmp_cell{sig_idx,int_idx};
                    else
                        collect_cell{j,k} = [collect_cell{j,k}, tmp_cell{sig_idx,int_idx}];
                    end
                    if isempty(recording_cell{j,k})
                        recording_cell{j,k} = i;
                    else
                        recording_cell{j,k} = [recording_cell{j,k} ,i];
                    end
                end                
            end
        end
    end
end

print_cell = {[]};
print_cell{size(collect_cell,1)+1,size(collect_cell,2)+1} = [];
print_cell(2:end,1) = signals';
print_cell{1,1} = 'number of values';
print_cell(1,2:end) = intervals;
for i = 1:size(collect_cell,1)
    for j = 1:size(collect_cell,2)
        print_cell{i+1,j+1} = length(collect_cell{i,j});
    end
end

writecell(print_cell,[app.settings.collect_path  app.edt_filename.Value '.xls' ],'FileType','spreadsheet','Sheet', 'general','Range','A1','WriteMode','replacefile');
app.lbl_working.Text = [num2str(round(100/(length(app.comb_stat)+1))) '% done'];

if ~isempty(app.comb_stat)
    for i = 1 : length(app.comb_stat)
        res = {[]};

        switch app.comb_stat(i).test
            case 'two sample t test'
                paired = true;
                test = 1;
            case 'Wilcoxon rank sum test'
                paired = true;
                test = 2;
            case 'anova'
                paired = false;
                test = 3;
            case 'anocova'
                paired = false;
                test = 4;
        end

        if paired
            int_idx(1) = find(contains(intervals,app.comb_stat(i).intervals{1}));
            int_idx(2) = find(contains(intervals,app.comb_stat(i).intervals{2}));
            res{1,1} = ['statistics group ' num2str(i)   ];
            res{1,2} = simple_name(app.comb_stat(i).intervals {1});
            res{1,3} = ' vs. ';
            res{1,4} = simple_name(app.comb_stat(i).intervals {2});
            res{2,1} = app.comb_stat(i).test;
            res{2,2} = 'null hypothesis';
            res{2,3} = 'p-value';
            
            for j = 1: length(app.comb_stat(i).signals)
                plot_pos = ['E' num2str((j-1)*2) '1'];
                sig_idx = find(contains(signals,app.comb_stat(i).signals{j}));
                switch test
                    case 1
                        [h,p] =  ttest2(collect_cell{sig_idx,int_idx(1)}, collect_cell{sig_idx,int_idx(2)});
                    case 2
                        [p,h] = ranksum(collect_cell{sig_idx,int_idx(1)}, collect_cell{sig_idx,int_idx(2)});
                end  
                res{end+1,1} = app.comb_stat(i).signals{j};
                res{end,2} = h;
                res{end,3} = p;
                stat_plot_cell = {collect_cell{sig_idx,int_idx(1)}, collect_cell{sig_idx,int_idx(2)}};
                combplot_boxplot_xls(stat_plot_cell, {res{1,2}, res{1,4}}, app.comb_stat(i).signals{j}, ...
                                 [app.settings.collect_path  app.edt_filename.Value '.xls' ],...
                                 ['statistics group ' num2str(i) ], h, p, plot_pos)
            end

        else
            
            names = [];
            namecell= {[]};

            for j = 1 : length(app.comb_stat(i).signals)
                
                if strcmp(app.comb_stat(i).posthoc,'none')
                    plot_pos = ['E' num2str((j-1)*2) '1'];
                else
                    plot_pos = ['H' num2str((j-1)*2) '1'];
                end

                vals = [];
                group = string([]);
                int_idx = app.comb_stat(i).intervals;
                tmp_cell = {[]};
                sig_idx = find(contains(signals,app.comb_stat(i).signals{j}));

                for k = 1 : length(app.comb_stat(i).intervals)
                    int_idx = find(contains(intervals,app.comb_stat(i).intervals{k}));

                    vals = [vals,collect_cell{sig_idx,int_idx}];
                    group(end+1:length(vals)) = string(simple_name(app.comb_stat(i).intervals{k}));
                    tmp_cell{k} = collect_cell{sig_idx,int_idx};
                    
                    if k == 1
                        names = [simple_name(app.comb_stat(i).intervals{k})];
                        varname = app.comb_stat(i).signals{j};
                    else
                        names = [names ' vs. ' simple_name(app.comb_stat(i).intervals{k})];
                    end
                    namecell{k} = simple_name(app.comb_stat(i).intervals{k});
                end

                switch app.comb_stat(i).test
                    case 'anova'
                        [p,~,stats] = anova1(vals,group,"off");
                    case 'anocova'
                        %[h,atab,ctab,stats] = aoctool(x,y,group,alpha,xname,yname,gname,"off")
                end
                                
                if strcmp(app.comb_stat(i).posthoc,'none')
                    if j  == 1
                        res{1,1} = ['statistics group ' num2str(i)   ];
                        res{1,2} = names; 
                        res{2,1} = app.comb_stat(i).test;
                        res{2,2} = 'p-value';
                    end
                    res{end+1,1} = varname;
                    res{end,2} = p;
                    combplot_boxplot_xls(tmp_cell,namecell , varname, ...
                                         [app.settings.collect_path  app.edt_filename.Value '.xls' ], ...
                                         ['statistics group ' num2str(i) ], false, p, plot_pos)
  
                else
                    c = multcompare(stats, 'CriticalValueType', app.comb_stat(i).posthoc);
                    p = c(c(:,6)<=0.05,[1,2,6]);
                    if j  == 1
                        res{1,1} = ['statistics group ' num2str(i)   ];
                        res{1,2} = names; 
                        res{1,3} = app.comb_stat(i).test;
                        res{1,4} = app.comb_stat(i).posthoc;
                        res{2,1} = 'int1';
                        res{2,2} = 'int2';
                        res{2,3} = 'mean difference';
                        res{2,4} = 'mean difference lower';
                        res{2,5} = 'mean difference upper';
                        res{2,6} = 'p-value';
                    end
                    res{end+1,1} = varname;
                    cc = num2cell(c);
                    for l = 1:size(c,1)
                        cc{l,1} = simple_name(app.comb_stat(i).intervals{c(l,1)});
                        cc{l,2} = simple_name(app.comb_stat(i).intervals{c(l,2)});
                    end
                    res(end+1:end+size(cc,1),1:6) = cc;
                    combplot_boxplot_phc_xls(tmp_cell,namecell , varname, ...
                                             [app.settings.collect_path  app.edt_filename.Value '.xls' ], ...
                                             ['statistics group ' num2str(i) ],  p, plot_pos)
                end
                
            end
        end

        writecell(res,[app.settings.collect_path  app.edt_filename.Value '.xls' ],'FileType','spreadsheet','Sheet', ['statistics group ' num2str(i) ],'Range','A1')
        app.lbl_working.Text = [num2str(round(100*(1+i)/(length(app.comb_stat)+1))) '% done'];
    end
end

end


% function output_str = simple_name(input_str)
% 
% tmp_idx = strfind(input_str, '<$§');
% if ~isempty(tmp_idx)
%     input_str(tmp_idx:end) = [];
% end
% tmp_idx = strfind(input_str,'§$>');
% if ~isempty(tmp_idx)
%     input_str(1:3) = [];
% end
% tmp_idx = strfind(input_str, '§$');
% if ~isempty(tmp_idx)
%     input_str(tmp_idx:tmp_idx+3) = [];
% end
% tmp_idx = strfind(input_str, '$§');
% if ~isempty(tmp_idx)
%     input_str(tmp_idx-2:tmp_idx+1) = [];
% end
% output_str = input_str;
% end

function combplot_boxplot_xls(datacell,namecell, ttl, filename,sheetname, hyp, p, startcell)
%%PLOT BOXPLOT AND insert into xls file

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

h = figure('Visible', "off");
boxplot2(datacell);
xticklabels(namecell)
title(ttl)
 
if hyp
    yl=ylim(gca);
    line([1,2],[yl(2)+.02*diff(yl) yl(2)+.02*diff(yl)], 'LineWidth', 2,'Color','k')
    ylim(gca,[yl(1) yl(2)+.1*diff(yl)])
    if p <= 0.001
        str = '***';
    elseif p <=0.01
        str = '**';
    elseif p <= 0.05
        str = '*';
    end
    text(1.5,yl(2)+.03*diff(yl),str,'FontSize',24,'HorizontalAlignment', 'center')
end

xlswritefig(h, filename, sheetname, startcell)
end

function combplot_boxplot_phc_xls(datacell,namecell, ttl, filename, sheetname,  p, startcell)
%%PLOT BOXPLOT AND insert into xls file

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

h = figure('Visible', "off");
boxplot2(datacell);
xticklabels(namecell)
title(ttl)
yl = ylim(gca);
dy = [.01*diff(yl),.02*diff(yl), .08*diff(yl)];

for i =1: size(p,1)
    line([p(i,1),p(i,2)],[yl(2)+dy(1), yl(2)+dy(1)], 'LineWidth', 2,'Color','k') 
    
    if p(i,3) <= 0.001
        str = '***';
    elseif p(i,3) <=0.01
        str = '**';
    elseif p(i,3) <= 0.05
        str = '*';
    end
    text(p(i,1)+diff(p(i,1:2))*0.5,yl(2)+dy(2),str,'FontSize',24,'HorizontalAlignment', 'center')
    yl(2) = yl(2)+dy(3);
end
ylim(gca,yl )
xlswritefig(h, filename, sheetname, startcell)
close(h)
end

