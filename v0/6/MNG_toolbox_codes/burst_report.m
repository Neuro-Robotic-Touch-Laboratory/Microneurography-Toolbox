
function burst_report(app)
%BURST_REPORT writes graphs and xls sheets with results of burst analysis
%   Detailed explanation goes here


%% create results

plot_data = struct('name',[],'integral',[],'amplitude',[],'duration',[],...
                   'n_bursts',[], 'burstrate',[], 'burstincidence',[], ...
                   'mean_hr',[]);

for i = 1 : length(app.burst_ints)
    plot_data(i).name = app.burst_ints(i).name;%plot_data(i).name = app.burst_ints(sorted_idx(i)).name;
    tmp_idx_burst = app.burst_res.use_burst(:,1) & app.burst_res.use_burst(:,2) & app.burst_res.use_burst(:,i+2);
    dur = sum(diff(app.burst_ints(i).borders,1,2));%dur = sum(diff(app.burst_ints(sorted_idx(i)).borders,1,2));
    plot_data(i).integral = app.burst_res.burst_int(tmp_idx_burst);
    plot_data(i).amplitude = app.burst_res.burst_amp(tmp_idx_burst);
    plot_data(i).duration = app.burst_res.burst_dur(tmp_idx_burst);
    plot_data(i).n_bursts = sum(tmp_idx_burst);
    plot_data(i).burstrate = sum(tmp_idx_burst)/dur;
    plot_data(i).dur = dur;
    
    if ~isempty(app.hb_res)
        tmp_idx_beats = app.hb_res.use_beats(:,1) & app.hb_res.use_beats(:,i+1);
        plot_data(i).burstincidence = sum(tmp_idx_burst)/sum(tmp_idx_beats);
        plot_data(i).mean_hr = sum(tmp_idx_beats)/dur *60;
        plot_data(i).rr_interval = app.hb_res.dt_instantaneous(tmp_idx_beats(1:end-1));
        plot_data(i).ts_hb = app.hb_res.t_events(tmp_idx_beats);
    end

end


%% plot/save boxplots

namecell = {[]};
tmp_cell_1 = {[]};
tmp_cell_2 = {[]};
tmp_cell_3 = {[]};
tmp_cell_4 = {[]};
tmp_cell_5 = {[]};
tmp_cell_6 = {[]};
tmp_cell_7 = {[]};
tmp_cell_8 = {[]};
for i = 1:length(plot_data)
    namecell{i} = plot_data(i).name;
    tmp_cell_1{i} = plot_data(i).integral;
    tmp_cell_2{i} = plot_data(i).amplitude;
    tmp_cell_3{i} = plot_data(i).duration;
    
    if ~isempty(app.hb_res)
        tmp_cell_4{i} = plot_data(i).rr_interval;
        tmp_cell_5{i} = HRV.HR(plot_data(i).rr_interval,10);
        tmp_cell_5{i}(isnan(tmp_cell_5{i})) = [];
        tmp_cell_6{i} = HRV.SDNN(plot_data(i).rr_interval,0);
        tmp_cell_7{i} = HRV.RMSSD(plot_data(i).rr_interval,0);
        tmp_cell_8{i} = HRV.pNN50(plot_data(i).rr_interval,0);
    end
end
%% get filename and replace '.' by '-'

path = app.settings.output_dir;
[indx,~] = listdlg('ListString',{'fig','jpg','eps'},'PromptString','select output file format');

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';
%%

plot_boxplot(tmp_cell_1,namecell,path,[file  '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_burst_analysis_burst_integral'],indx,'Burst Integral')

plot_boxplot(tmp_cell_2,namecell,path,[file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_burst_analysis_burst_amplitude'],indx,'Burst Amplitude')

plot_boxplot(tmp_cell_3,namecell,path,[file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_burst_analysis_burst_duration'],indx, 'Burst Duration')
if ~isempty(app.hb_res)
    plot_boxplot(tmp_cell_4,namecell,path,[file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_burst_analysis_rr-interval'],indx, 'RR-Interval')
end


%% calculate / organize results to be written to xls file

if ~isempty(app.bp_res)
    sys_ch_idx = find(contains(vertcat(app.data.name),'systolic BP'));
    mea_ch_idx = find(contains(vertcat(app.data.name),'mean BP'));
    dia_ch_idx = find(contains(vertcat(app.data.name),'diastolic BP'));
else
    sys_ch_idx = nan;
    mea_ch_idx = nan;
    dia_ch_idx = nan;
end

if ~isempty(app.resp_res)
    resp_ch_idx = find(contains(vertcat(app.data.name),'Respiration Rate'));
else
    resp_ch_idx = nan;
end
if ~isempty(find(contains(vertcat(app.data.name),'etCO2')))
    co2_ch_idx = find(contains(vertcat(app.data.name),'etCO2'));
else
    co2_ch_idx = nan;
end

plot_cell = {[]};
plot_cell{1,1} = 'interval';
plot_cell{2,1} = 'duration [s]';
plot_cell{3,1} = 'number of bursts';
plot_cell{4,1} = 'burst rate [Hz]';
plot_cell{5,1} = 'median integral [uV*ms]';
plot_cell{6,1} = 'mean integral [uV*ms]';
plot_cell{7,1} = 'std integral [uV*ms]';
plot_cell{8,1} = 'median amplidude [uV]';
plot_cell{9,1} = 'mean amplidude [uV]';
plot_cell{10,1} = 'std amplidude [uV]';
plot_cell{11,1} = 'median duration [s]';
plot_cell{12,1} = 'mean duration [s]';
plot_cell{13,1} = 'std duration [s]';

if ~isempty(app.hb_res)
    hb_idx = size(plot_cell,1)+1;
    plot_cell{hb_idx,1} = 'burst incidence';
    plot_cell{hb_idx+1,1} = 'mean heart rate [BPM]';
    plot_cell{hb_idx+2,1} = 'median rr-interval [s]';
    plot_cell{hb_idx+3,1} = 'mean rr-interval [s]';
    plot_cell{hb_idx+4,1} = 'std rr-interval [s]';
    plot_cell{hb_idx+5,1} = 'median heartrate [BPM]';
    plot_cell{hb_idx+6,1} = 'mean heartrate [BPM]';
    plot_cell{hb_idx+7,1} = 'std heartrate [BPM]';
    plot_cell{hb_idx+8,1} = 'SDNN [ms]';
    plot_cell{hb_idx+9,1} = 'RMSSD [ms]';
    plot_cell{hb_idx+10,1} = 'pNN50 [%]';
    plot_cell{hb_idx+11,1} = 'pLF [%]';
    plot_cell{hb_idx+12,1} = 'pHF [%]';
    plot_cell{hb_idx+13,1} = 'LF/HF ratio';
    plot_cell{hb_idx+14,1} = 'VLF [ms²]';
    plot_cell{hb_idx+15,1} = 'LF [ms²]';
    plot_cell{hb_idx+16,1} = 'HF [ms²]';
end

if ~isempty(app.bp_res)
    bp_idx = size(plot_cell,1)+1;
    plot_cell{bp_idx,1} = 'median systolic BP [mmHg]';
    plot_cell{bp_idx+1,1} = 'mean systolic BP [mmHg]';
    plot_cell{bp_idx+2,1} = 'std systolic BP [mmHg]';
    plot_cell{bp_idx+3,1} = 'median mean BP [mmHg]';
    plot_cell{bp_idx+4,1} = 'mean mean BP [mmHg]';
    plot_cell{bp_idx+5,1} = 'std mean BP [mmHg]';
    plot_cell{bp_idx+6,1} = 'median diastolic BP [mmHg]';
    plot_cell{bp_idx+7,1} = 'mean diastolic BP [mmHg]';
    plot_cell{bp_idx+8,1} = 'std diastolic BP [mmHg]';  
else
    bp_idx = nan;
end

if ~isempty(app.resp_res)
    resp_idx = size(plot_cell,1)+1;
    plot_cell{resp_idx,1} = 'median respiration rate [BPM]';
    plot_cell{resp_idx+1,1} = 'mean respiration rate [BPM]';
    plot_cell{resp_idx+2,1} = 'std respiration rate [BPM]';
else
    resp_idx = nan;
end
%% baroreflexstuff
if ~isempty(app.bp_res) && ~isempty(app.hb_res)
    baro_idx = size(plot_cell,1)+1;
    plot_cell{baro_idx} = 'std baroreflexsensitivity [ms/mmHg]';
else
    baro_idx = nan;
end

if ~isnan(co2_ch_idx)
    etco2_idx = size(plot_cell,1)+1;
    plot_cell{etco2_idx,1} = 'median etCO2 [mmHg]';
    plot_cell{etco2_idx+1,1} = 'mean etCO2 [mmHg]';
    plot_cell{etco2_idx+2,1} = 'std etCO2 [mmHg]';
else
    etco2_idx = nan;
end

if ~isempty(app.resp_res)
    minvent_idx = size(plot_cell,1)+1;
    plot_cell{resp_idx,1} = 'median minute ventilation';
    plot_cell{resp_idx+1,1} = 'mean minute ventilation';
    plot_cell{resp_idx+2,1} = 'std minute ventilation';
else
    minvent_idx = nan;
end

for i = 1: length(plot_data)
    plot_cell{1,i+1} = plot_data(i).name;
    plot_cell{2,i+1} = plot_data(i).dur;
    plot_cell{3,i+1} = plot_data(i).n_bursts;
    plot_cell{4,i+1} = plot_data(i).burstrate;
    plot_cell{5,i+1} = median(plot_data(i).integral);
    plot_cell{6,i+1} = mean(plot_data(i).integral);
    plot_cell{7,i+1} = std(plot_data(i).integral);
    plot_cell{8,i+1} = median(plot_data(i).amplitude);
    plot_cell{9,i+1} = mean(plot_data(i).amplitude);
    plot_cell{10,i+1} = std(plot_data(i).amplitude);
    plot_cell{11,i+1} = median(plot_data(i).duration);
    plot_cell{12,i+1} = mean(plot_data(i).duration);
    plot_cell{13,i+1} = std(plot_data(i).duration);
   
    if ~isempty(app.hb_res)
        plot_cell{hb_idx,i+1} = plot_data(i).burstincidence;
        plot_cell{hb_idx+1,i+1} = plot_data(i).mean_hr;
        plot_cell{hb_idx+2,i+1} = median(plot_data(i).rr_interval);
        plot_cell{hb_idx+3,i+1} = mean(plot_data(i).rr_interval);
        plot_cell{hb_idx+4,i+1} = std(plot_data(i).rr_interval);
        plot_cell{hb_idx+5,i+1} = median(tmp_cell_5{i});
        plot_cell{hb_idx+6,i+1} = mean(tmp_cell_5{i});
        plot_cell{hb_idx+7,i+1} = std(tmp_cell_5{i});
        plot_cell{hb_idx+8,i+1} = tmp_cell_6{i};
        plot_cell{hb_idx+9,i+1} = tmp_cell_7{i};
        plot_cell{hb_idx+10,i+1} = tmp_cell_8{i};
        [pLF,pHF,LFHFratio,VLF,LF,HF,~,~,~] = ...
            HRV.fft_val_fun(plot_data(i).rr_interval, 1/app.data(app.settings.channel_idx.ecg).ts(1),'spline');
        plot_cell{hb_idx+11,i+1} = pLF;
        plot_cell{hb_idx+12,i+1} = pHF;
        plot_cell{hb_idx+13,i+1} = LFHFratio;
        plot_cell{hb_idx+14,i+1} = VLF;
        plot_cell{hb_idx+15,i+1} = LF;
        plot_cell{hb_idx+16,i+1} = HF;
    end
    if ~isnan(bp_idx)
        tmp_sig = [];
        tmp_ts = app.data(sys_ch_idx).ts(1);
        for j = 1:size(app.burst_ints(i).borders,1)
            tmp_sig = [tmp_sig, [app.data(sys_ch_idx).data(app.burst_ints(i).borders(j,1)/tmp_ts:app.burst_ints(i).borders(j,2)/tmp_ts)';
                                 app.data(mea_ch_idx).data(app.burst_ints(i).borders(j,1)/tmp_ts:app.burst_ints(i).borders(j,2)/tmp_ts)';
                                 app.data(dia_ch_idx).data(app.burst_ints(i).borders(j,1)/tmp_ts:app.burst_ints(i).borders(j,2)/tmp_ts)']];
        end
        plot_cell{bp_idx,i+1} = median(tmp_sig(1,:));
        plot_cell{bp_idx+1,i+1} = mean(tmp_sig(1,:));
        plot_cell{bp_idx+2,i+1} = std(tmp_sig(1,:));
        plot_cell{bp_idx+3,i+1} = median(tmp_sig(2,:));
        plot_cell{bp_idx+4,i+1} = mean(tmp_sig(2,:));
        plot_cell{bp_idx+5,i+1} = std(tmp_sig(2,:));
        plot_cell{bp_idx+6,i+1} = median(tmp_sig(3,:));
        plot_cell{bp_idx+7,i+1} = mean(tmp_sig(3,:));
        plot_cell{bp_idx+8,i+1} = std(tmp_sig(3,:));  
    end

    if ~isnan(resp_idx)
        tmp_sig = [];
        tmp_ts = app.data(resp_ch_idx).ts(1);
        for j = 1:size(app.burst_ints(i).borders,1)
            tmp_sig = [tmp_sig, app.data(resp_ch_idx).data(app.burst_ints(i).borders(j,1)/tmp_ts:app.burst_ints(i).borders(j,2)/tmp_ts)'];
        end
        plot_cell{resp_idx,i+1} = median(tmp_sig);
        plot_cell{resp_idx+1,i+1} = mean(tmp_sig);
        plot_cell{resp_idx+2,i+1} = std(tmp_sig); 
    end

    if ~isnan(baro_idx)
        plot_cell{baro_idx,i+1} = plot_cell{hb_idx+4,i+1}*1000 / plot_cell{bp_idx+2,i+1};
    end


    if ~isnan(etco2_idx)
        tmp_sig = [];
        tmp_ts = app.data(co2_ch_idx).ts(1);
        for j = 1:size(app.burst_ints(i).borders,1)
            tmp_sig = [tmp_sig, app.data(co2_ch_idx).data(app.burst_ints(i).borders(j,1)/tmp_ts:app.burst_ints(i).borders(j,2)/tmp_ts)'];
        end
        plot_cell{etco2_idx,i+1} = median(tmp_sig);
        plot_cell{etco2_idx+1,i+1} = mean(tmp_sig);
        plot_cell{etco2_idx+2,i+1} = std(tmp_sig); 
    end
%     if ~isnan(minvent_idx)
%     
%     end

end

writecell(plot_cell,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],'FileType','spreadsheet','Sheet', 'general','Range','A1','WriteMode','replacefile');

app.lbl_working.Text = [num2str(round((1) / (1+length(app.stat_group))*100)) '% done'];
pause(.05)

if ~isempty(app.stat_group)
    for i = 1 : length(app.stat_group)
        res = {[]};
        switch app.stat_group(i).paired
            case true
                switch app.stat_group(i).test
                    case 'two sample t test'
                        [h,p] =  ttest2(plot_data(app.stat_group(i).intervals(1)).integral, plot_data(app.stat_group(i).intervals(2)).integral);
                        test_string = 'two sample t test';
                    case 'Wilcoxon rank sum test'
                        [p,h] = ranksum(plot_data(app.stat_group(i).intervals(1)).integral, plot_data(app.stat_group(i).intervals(2)).integral);
                        test_string = 'Wilcoxon rank sum test';
                end
        
                res{1,1} = ['statistics group ' num2str(i)   ];
                res{1,2} = plot_data(app.stat_group(i).intervals(1)).name; 
                res{1,3} = ' vs. ';
                res{1,4} = plot_data(app.stat_group(i).intervals(2)).name; 
                res{2,1} = test_string;
                res{2,2} = 'null hypothesis';
                res{2,3} = 'p-value';
                res{3,1} = 'Burst Integral';
                if islogical(h)
                    h = int8(h);
                end
                res{3,2} = h;
                res{3,3} = p;
                stat_plot_cell = {plot_data(app.stat_group(i).intervals(1)).integral, plot_data(app.stat_group(i).intervals(2)).integral};
                plot_boxplot_xls(stat_plot_cell, {plot_data(app.stat_group(i).intervals(1)).name, plot_data(app.stat_group(i).intervals(2)).name}, ...
                                 'Burst Integral', [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ], ['statistics group ' num2str(i) ], h, p, 'E3')
                
                switch app.stat_group(i).test
                    case 'two sample t test'
                        [h,p] =  ttest2(plot_data(app.stat_group(i).intervals(1)).amplitude, plot_data(app.stat_group(i).intervals(2)).amplitude);
                    case 'Wilcoxon rank sum test'
                        [p,h] = ranksum(plot_data(app.stat_group(i).intervals(1)).amplitude, plot_data(app.stat_group(i).intervals(2)).amplitude);
                end
                res{4,1} = 'Burst Amplitude';
                if islogical(h)
                    h = int8(h);
                end
                res{4,2} = h;
                res{4,3} = p;
                stat_plot_cell = {plot_data(app.stat_group(i).intervals(1)).amplitude, plot_data(app.stat_group(i).intervals(2)).amplitude};
                plot_boxplot_xls(stat_plot_cell, {plot_data(app.stat_group(i).intervals(1)).name, plot_data(app.stat_group(i).intervals(2)).name},... 
                                 'Burst Amplitude', [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ], ['statistics group ' num2str(i) ], h, p, 'k3')
                switch app.stat_group(i).test
                    case 'two sample t test'
                        [h,p] =  ttest2(plot_data(app.stat_group(i).intervals(1)).duration, plot_data(app.stat_group(i).intervals(2)).duration);
                    case 'Wilcoxon rank sum test'
                        [p,h] = ranksum(plot_data(app.stat_group(i).intervals(1)).duration, plot_data(app.stat_group(i).intervals(2)).duration);
                end
                res{5,1} = 'Burst Duration';
                if islogical(h)
                    h = int8(h);
                end
                res{5,2} = h;
                res{5,3} = p; 
                stat_plot_cell = {plot_data(app.stat_group(i).intervals(1)).duration, plot_data(app.stat_group(i).intervals(2)).duration};
                plot_boxplot_xls(stat_plot_cell, {plot_data(app.stat_group(i).intervals(1)).name, plot_data(app.stat_group(i).intervals(2)).name}, 'Burst Duration', [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],...
                                 ['statistics group ' num2str(i) ], h, p, 'Q3')
                if ~isempty(app.hb_res)
                    switch app.stat_group(i).test
                        case 'two sample t test'
                            [h,p] =  ttest2(tmp_cell_4{app.stat_group(i).intervals(1)}, tmp_cell_4{app.stat_group(i).intervals(2)});
                        case 'Wilcoxon rank sum test'
                            [p,h] = ranksum(tmp_cell_4{app.stat_group(i).intervals(1)}, tmp_cell_4{app.stat_group(i).intervals(2)});
                    end
                    res{6,1} = 'RR-intervals';
                    if islogical(h)
                        h = int8(h);
                    end
                    res{6,2} = h;
                    res{6,3} = p; 
                    stat_plot_cell = {tmp_cell_4{app.stat_group(i).intervals(1)}, tmp_cell_4{app.stat_group(i).intervals(2)}};
                    plot_boxplot_xls(stat_plot_cell, {plot_data(app.stat_group(i).intervals(1)).name, plot_data(app.stat_group(i).intervals(2)).name}, 'RR-intervals', [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],...
                                     ['statistics group ' num2str(i) ], h, p, 'W3')
                    
                    switch app.stat_group(i).test
                        case 'two sample t test'
                            [h,p] =  ttest2(tmp_cell_5{app.stat_group(i).intervals(1)}, tmp_cell_5{app.stat_group(i).intervals(2)});
                        case 'Wilcoxon rank sum test'
                            [p,h] = ranksum(tmp_cell_5{app.stat_group(i).intervals(1)}, tmp_cell_5{app.stat_group(i).intervals(2)});
                    end
                    res{7,1} = 'averaged heart rate';
                    if islogical(h)
                        h = int8(h);
                    end
                    res{7,2} = h;
                    res{7,3} = p;
                    stat_plot_cell = {tmp_cell_5{app.stat_group(i).intervals(1)}, tmp_cell_5{app.stat_group(i).intervals(2)}};
                    plot_boxplot_xls(stat_plot_cell, {plot_data(app.stat_group(i).intervals(1)).name, plot_data(app.stat_group(i).intervals(2)).name}, 'averaged heart rate', [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],...
                                     ['statistics group ' num2str(i) ], h, p, 'AD3')
                end
                writecell(res,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],'FileType','spreadsheet','Sheet', ['statistics group ' num2str(i) ],'Range','A1')
    
            case false
                 names = [];
                 namecell= {[]};
                 if ~isempty(app.hb_res)
                     signals = 1:5;
                 else
                     signals = 1:3;
                 end
                 for j = signals
                    vals = [];
                    group = string([]);
                    
                    int_idx = app.stat_group(i).intervals;
                    tmp_cell = {[]};
                    for k = 1:length(int_idx)
                        switch j
                            case 1
                                vals = [vals,plot_data(int_idx(k)).integral];
                                group(end+1:length(vals)) = string(app.burst_ints(int_idx(k)).name);
                                tmp_cell{k} = plot_data(int_idx(k)).integral;
                                if k == 1
                                    names = [app.burst_ints(int_idx(k)).name];
                                    varname = 'Burst integral';
                                    startcell = {'E1','H1'};
                                else
                                    names = [names ' vs. ' app.burst_ints(int_idx(k)).name];
                                end
                                namecell{k} = app.burst_ints(int_idx(k)).name;
                                
                            case 2
                                vals = [vals,plot_data(int_idx(k)).amplitude];
                                group(end+1:length(vals)) = string(app.burst_ints(int_idx(k)).name);
                                tmp_cell{k} = plot_data(int_idx(k)).amplitude;
                                if k == 1
                                    varname = 'Burst amplitude';
                                    startcell = {'E21', 'H23'};
                                end
                                
                            case 3
                                vals = [vals,plot_data(int_idx(k)).duration'];
                                group(end+1:length(vals)) = string(app.burst_ints(int_idx(k)).name);
                                tmp_cell{k} = plot_data(int_idx(k)).duration;
                                if k == 1
                                    varname = 'Burst duration';
                                    startcell = {'E41','H45'};
                                end
                            case 4
                                vals = [vals,tmp_cell_4{int_idx(k)}'];
                                group(end+1:length(vals)) = string(app.burst_ints(int_idx(k)).name);
                                tmp_cell{k} = tmp_cell_4{int_idx(k)};
                                if k == 1
                                    varname = 'RR-intervals';
                                    startcell = {'E61','H67'};
                                end
                                
                            case 5
                                vals = [vals,tmp_cell_5{int_idx(k)}'];
                                group(end+1:length(vals)) = string(app.burst_ints(int_idx(k)).name);
                                tmp_cell{k} = tmp_cell_5{int_idx(k)};
                                if k == 1
                                    varname = 'averaged heart rate';
                                    startcell = {'E81', 'H89'};
                                end     
                        end
                    end
                    
                    switch app.stat_group(i).test
                        case 'anova'
                            [p,~,stats] = anova1(vals,group,"off");
                        case 'anocova'
                            %[h,atab,ctab,stats] = aoctool(x,y,group,alpha,xname,yname,gname,"off")
                    end
                                    
                    if strcmp(app.stat_group(i).posthoc , 'none')
                        if j  == 1
                            res{1,1} = ['statistics group ' num2str(i)   ];
                            res{1,2} = names; 
                            res{2,1} = app.stat_group(i).test;
                            res{2,2} = 'p-value';
                        end
                        res{2+j,1} = varname;
                        res{2+j,2} = p;
                        plot_boxplot_xls(tmp_cell,namecell , varname, [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ], ['statistics group ' num2str(i) ], false, p, startcell{1})
      
                    else
                        c = multcompare(stats, 'CriticalValueType', app.stat_group(i).posthoc);
                        p = c(c(:,6)<=0.05,[1,2,6]);
                        if j  == 1
                            res{1,1} = ['statistics group ' num2str(i)   ];
                            res{1,2} = names; 
                            res{1,3} = app.stat_group(i).test;
                            res{1,4} = app.stat_group(i).posthoc;
                            res{2,1} = 'int1';
                            res{2,2} = 'int2';
                            res{2,3} = 'mean difference';
                            res{2,4} = 'mean difference lower';
                            res{2,5} = 'mean difference upper';
                            res{2,6} = 'p-value';
                        end
                        res{1+((size(c,1)+1)*(j-1))+2,1} = varname;
                        cc = num2cell(c);
                        for k = 1:size(c,1)
                            cc{k,1} = app.burst_ints(int_idx(c(k,1))).name;
                            cc{k,2} = app.burst_ints(int_idx(c(k,2))).name;
                        end
                        res(2+((size(c,1)+1)*(j-1))+2: 1+((size(c,1)+1)*(j-1))+2+size(cc,1),1:6) = cc; 
                        plot_boxplot_phc_xls(tmp_cell,namecell , varname, [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ], ['statistics group ' num2str(i) ],  p, startcell{2})
                        
                    end
    
                end
                writecell(res,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],'FileType','spreadsheet','Sheet', ['statistics group ' num2str(i) ],'Range','A1')
            
        end
        app.lbl_working.Text = [num2str(round((1+i) / (1+length(app.stat_group))*100)) '% done'];
        pause(.05)
    end
end


end


function plot_boxplot(datacell,namecell,path,file,indx, ttl)
%%PLOT BOXPLOT AND SAVE TO FILE 

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

h = figure('Visible', "off");
boxplot2(datacell);
xticklabels(namecell)
title(ttl)
for i = 1 : length(indx)
    switch indx(i)
        case 1
            %h.Visible = 'on';
            savefig(h,[path '\' file '.fig'],'compact')
            %switch_vis([path '\' file '.fig'])
        case 2
            saveas(h,[path '\' file '.jpeg'])
        case 3
            saveas(h,[path '\' file '.epsc'])
    end
end
close (h)
end

function plot_boxplot_xls(datacell,namecell, ttl, filename,sheetname, hyp, p, startcell)
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

function plot_boxplot_phc_xls(datacell,namecell, ttl, filename, sheetname,  p, startcell)
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
