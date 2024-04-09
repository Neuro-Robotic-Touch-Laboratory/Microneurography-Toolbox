function calc_combine_stats(app)

intervals = {[]};
signals = {[]};
anocovas= [];

for i = 1: length(app.comb_stat)
    if i == 1
        intervals = app.comb_stat(i).intervals;
        signals = app.comb_stat(i).signals;
    else
        intervals = unique(horzcat(intervals, app.comb_stat(i).intervals),'stable');
        signals = unique(horzcat(signals, app.comb_stat(i).signals),'stable');
    end
    if strcmp(app.comb_stat(i).test, 'anocova')
        anocovas = [anocovas,i];
    end
end
use_rec = app.settings.brst_files(contains(app.lstbx_rec_combine.Items,app.lstbx_rec_combine.Value));

collect_cell = {[]};
collect_cell{length(signals),length(intervals)} = [];
l = length(use_rec);
for i=1:length(signals)
    for j=1:length(intervals)
        collect_cell{i,j} = nan(1,l);
    end
end



for i = 1: length(use_rec)
    tmp_cell = readcell(use_rec(i).filename,'Sheet', 'general');
    for j = 1 : length(signals)
        tmp_idx = find(contains(use_rec(i).signals(:,1),signals{j}));
        if ~isempty(tmp_idx)
            sig_idx = app.settings.brst_files(i).signals{tmp_idx,2};
            
            for k = 1 : length(intervals)
                tmp_idx = find(contains(app.settings.brst_files(i).subintervals(:,1),intervals{k}));
                if ~isempty(tmp_idx)
                    int_idx = app.settings.brst_files(i).subintervals{tmp_idx,2};
                    collect_cell{j,k}(i) = tmp_cell{sig_idx,int_idx};
                end                
            end
        end
    end
end

covar = struct('name', [], 'values', [], 'asprev', nan,'checkarray',[], 'test_idx',[]);
rem_ano = [];

idx = 1;
for i = 1 : length(anocovas)
    tmp_chk = false(l,length(app.comb_stat(anocovas(i)).intervals));
    sig_idx = find(contains(signals,app.comb_stat(anocovas(i)).signals));
    int_idx = find(contains(intervals,app.comb_stat(anocovas(i)).intervals));
    for j = 1:length(sig_idx)
        for k = 1:length(int_idx)
            tmp_chk(:,k) = tmp_chk(:,k) | ~isnan((collect_cell{sig_idx(j),int_idx(k)})'); 
        end 
    end

    if ~isempty(find(sum(tmp_chk,2) >1))
        app.comb_stat(anocovas(i)).test = 'anova';
        anocovas(i) = nan;
    else
%         covar(idx).anoco = true;
        for j = 2:size(tmp_chk,2)
            tmp_chk(:,1) = or(tmp_chk(:,1),tmp_chk(:,j));
        end
        covar(idx).checkarray = tmp_chk(:,1);
        app.comb_stat(anocovas(i)).co_idx = idx;
        covar(idx).test_idx = anocovas(i);
        if idx > 1
            for j = idx-1:-1:1
                if isequal(covar(idx).checkarray,covar(j).checkarray)
                    covar(idx).asprev = j;
                end
            end
        end
        idx = idx+1;
    end
end

anocovas(isnan(anocovas)) = [];
if ~isempty(anocovas)
    covar = covardlg(covar,use_rec,anocovas);
end
%% Finish !!!

print_cell = {[]};
print_cell{size(collect_cell,1)+1,size(collect_cell,2)+1} = [];
print_cell(2:end,1) = signals';
print_cell{1,1} = 'number of values';
print_cell(1,2:end) = intervals;
for i = 1:size(collect_cell,1)
    for j = 1:size(collect_cell,2)
        print_cell{i+1,j+1} = length(~isnan(collect_cell{i,j}));
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
                val1 = collect_cell{sig_idx,int_idx(1)};
                val1(isnan(val1)) = [];
                val2 = collect_cell{sig_idx,int_idx(2)};
                val2(isnan(val2)) = [];
                switch test
                    case 1
                        [h,p] =  ttest2(val1, val2);
                    case 2
                        [p,h] = ranksum(val1, val2);
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
                cov = [];
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
                    if strcmp(app.comb_stat(i).test,'anocova')
                        cov = [cov;covar(app.comb_stat(i).co_idx).values];
                    end
                end
                group(isnan(vals)) = [];
                if strcmp(app.comb_stat(i).test,'anocova')
                    cov(isnan(vals)) = [];
                end
                vals(isnan(vals)) = [];
                switch app.comb_stat(i).test
                    case 'anova'
                        [p,~,stats] = anova1(vals,group,"off");
                    case 'anocova'
                        [h,atab,ctab,stats] = aoctool(vals,cov,group,0.05,'','','',"off");
                end
                error_flag = false;
                if ~strcmp(app.comb_stat(i).posthoc,'none')
                    try
                        c = multcompare(stats, 'CriticalValueType', app.comb_stat(i).posthoc, 'Display','off');
                    catch ME
                        error_flag = true;
                        error_msg = [ME.message ' - Posthoc failed'];
                    end
                end

                if strcmp(app.comb_stat(i).posthoc,'none') || error_flag 
                    if j  == 1
                        res{1,1} = ['statistics group ' num2str(i)   ];
                        res{1,2} = names; 
                        res{2,1} = app.comb_stat(i).test;
                        res{2,2} = 'p-value';
                    end
                    res{end+1,1} = varname;
                    
                    if error_flag
                        res{end,2} = error_msg;
                    else
                        res{end,2} = p;
                    end
                    combplot_boxplot_xls(tmp_cell,namecell , varname, ...
                                         [app.settings.collect_path  app.edt_filename.Value '.xls' ], ...
                                         ['statistics group ' num2str(i) ], false, p, plot_pos)
  
                else
                    
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

