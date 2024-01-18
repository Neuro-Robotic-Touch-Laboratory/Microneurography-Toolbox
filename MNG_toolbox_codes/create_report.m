function create_report(app)

tmp = sort(app.lstbx_report.Value);
idx = nan(length(tmp),2);
for i = 1: length(tmp)
    tmp_idx = [find(tmp{i} =='[',1),find(tmp{i} =='-',1), find(tmp{i} ==']',1)];
    idx(i,1) = int16(str2num(tmp{i}(tmp_idx(1)+1:tmp_idx(2)-1)));
    idx(i,2) = int16(str2num(tmp{i}(tmp_idx(2)+1:tmp_idx(3)-1)));
end
tbl_figs ={[]};
filenames = {[]};
cell_idx = 1;

tmp = unique (idx(:,1));
filename = 'report';
for i = 1:length(tmp)
    filename = [filename '_' app.settings.figs(tmp(i)).name '.pdf'];
end
filename = [app.settings.rep_dir '\' filename];

if isfile(filename)
    delete(filename)
end

for rec_int = 1: size(idx,1)

    rec_int_name =  [app.settings.figs(idx(rec_int,1)).name '_' app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name];
    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.burst.done
        tbl_figs{cell_idx,1} = [char(rec_int_name) ' burstanalysis - burst amplitude overview'];
        filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                 char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                 char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name)...
                                 '_burst_analysis_burst_amplitude.'];
        filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.burst.amp.format));
        cell_idx = cell_idx+1;
        
        tbl_figs{cell_idx,1} = [char(rec_int_name) ' burstanalysis - burst duration overview'];
        filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                 char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                 char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name)...
                                 '_burst_analysis_burst_duration.'];
        filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.burst.dur.format));
        cell_idx = cell_idx+1;
        
        tbl_figs{cell_idx,1} = [char(rec_int_name) ' burstanalysis - burst integral overview'];
        filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                 char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                 char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name)...
                                 '_burst_analysis_burst_integral.'];
        filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.burst.int.format));
        cell_idx = cell_idx+1;

        tbl_figs{cell_idx,1} = [char(rec_int_name) ' burstanalysis - rr-interval overview'];
        filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                 char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                 char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name)...
                                 '_burst_analysis_rr-interval.'];
        filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.burst.rr.format));
        cell_idx = cell_idx+1;
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint)
            tbl_figs{cell_idx,1} = [char(rec_int_name) ' spikeanalysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).name) ' overview'];
            filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                     char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_spike_analysis_overview_'...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).name) '.'];
            filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).format));
            cell_idx = cell_idx+1;

            for j = 1:length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).unit)
                tbl_figs{cell_idx,1} = [char(rec_int_name) ' spikeanalysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).name) ' unit ' num2str(j)];
                filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                         char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_spike_analysis_unit_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).unit(j).name) '_'...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).name) '.'];
                filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spike.subint(i).unit(j).format));
                cell_idx = cell_idx+1;
            end
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint)
            for j = 1:length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).signal)
                tbl_figs{cell_idx,1} = [char(rec_int_name) ' spectralanalysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).name) ' - ' ...
                                      char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).signal(j).name)];
                filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                         char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_spectral_analysis_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).name) '_SIG_'...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).signal(j).name) '.'];
                filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.spectral.subint(i).signal(j).format));
                cell_idx = cell_idx+1;
            end
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint)
            for j = 1:length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).signal)
                tbl_figs{cell_idx,1} = [char(rec_int_name) ' waveletanalysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).name) ' - ' ...
                                      char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).signal(j).name)];
                filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                         char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_wavelet_analysis_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).name) '_SIG_'...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).signal(j).name) '.'];
                filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.wavelet.subint(i).signal(j).format));
                cell_idx = cell_idx+1;
            end
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint)
            for j = 1:length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).signal)
                tbl_figs{cell_idx,1} = [char(rec_int_name) ' entropy analysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).name) ' - ' ...
                                      char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).signal(j).name)];
                filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                         char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_entropy_analysis_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).name) '_SIG_'...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).signal(j).name) '.'];
                filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.entropy.subint(i).signal(j).format));
                cell_idx = cell_idx+1;
            end
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint)
            for j = 1:length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).signal)
                tbl_figs{cell_idx,1} = [char(rec_int_name) ' correlationanalysis - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).name) ' - ' ...
                                      char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).signal(j).name)];
                filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                         char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_correlation_' ...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).name) '_SIG_'...
                                         char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).signal(j).name) '.'];
                filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.correlation.subint(i).signal(j).format));
                cell_idx = cell_idx+1;
            end
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.stepwise.done
         for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.stepwise.subint)
            tbl_figs{cell_idx,1} = [char(rec_int_name) ' stepwise regression - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.stepwise.subint(i).name)];
            filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                     char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_stepwise_regression_'...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.stepwise.subint(i).name) '.'];
            filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.stepwise.subint(i).format));
            cell_idx = cell_idx+1;
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.annotate.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.annotate.subint)
            tbl_figs{cell_idx,1} = [char(rec_int_name) ' annotate - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.annotate.subint(i).name)];
            filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                     char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_annotate_'...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.annotate.subint(i).name) '.'];
            filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.annotate.subint(i).format));
            cell_idx = cell_idx+1;
        end
    end

    if app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.transduction.done
        for i =1: length(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.transduction.subint)
            tbl_figs{cell_idx,1} = [char(rec_int_name) ' transduction - ' char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.transduction.subint(i).name)];
            filenames{cell_idx,1} = [char(app.settings.rep_dir) '\' ...
                                     char(app.settings.figs(idx(rec_int,1)).name) '_INT_' ...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).name) '_transduction_'...
                                     char(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.transduction.subint(i).name) '.'];
            filenames{cell_idx,2} = min(find(app.settings.figs(idx(rec_int,1)).ints(idx(rec_int,2)).anas.transduction.subint(i).format));
            cell_idx = cell_idx+1;
        end
    end
end

num_figs = length(tbl_figs);
dp = ceil(num_figs/70);
for i= 1:num_figs
    tbl_figs{i,2} = i+dp;
end

while ~isempty(tbl_figs)
    tbl_data = {[]};

    if size(tbl_figs,1) >=40
        l =40;
    else 
        l = size(tbl_figs,1);
    end
    tbl_data(1:l,1:2) = tbl_figs(1:l,:);
    tbl_figs(1:l,:) = [];

    if size(tbl_figs,1) >=40
        l =40;
    else 
        l = size(tbl_figs,1);
    end
    tbl_data(1:l,3:4) = tbl_figs(1:l,:);
    tbl_figs(1:l,:) = [];

    h = uifigure('Position', get(0, 'Screensize'),'Visible', 'off');
    tbl = uitable(h,'Data',tbl_data,'ColumnName',{'Figure','Page', 'Figure', 'Page'},'RowName',[], 'Position', [0,0,h.Position(3),h.Position(4)]);
    if ~isfile(filename)
        exportapp(h,filename)
    else
        exportapp(h,'tmp.pdf')
        mergePdfs({filename, 'tmp.pdf'}, filename)
    end
    close(h)
end

for i= 1:size(filenames,1)
    switch filenames{i,2}
        case 1
            h =  openfig([filenames{i,1} 'fig'],'invisible');
            exportgraphics(h,filename,'Append',true)
            close(h)
        case 2
            im = imread([filenames{i,1} 'jpeg']);
            h = figure('Position', get(0, 'Screensize'),'Visible', 'off');
            imshow(im, 'Parent',gca);
            exportgraphics(h,filename,'Append',true)
            close(h)
        case 3
            h = figure('Visible', 'off');
            text(gca,.3,.5, 'eps files cannot be loaded', 'FontWeight', 'bold')
            exportgraphics(h,filename,'Append',true)
            close(h)
    end
    app.lbl_working.Text = [num2str(i) ' of ' num2str(num_figs) ' appended' ];
end

   
end