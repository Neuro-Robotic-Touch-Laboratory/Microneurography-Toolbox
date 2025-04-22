
function burst_report(app)
%BURST_REPORT writes graphs and xls sheets with results of burst analysis
%   Detailed explanation goes here
tmp = [];
for i = 1 :length(app.stat_group)
    tmp = [tmp ,app.stat_group(i).signals];
end
tmp = unique(tmp);
name_str = {};
for i = 1 : length(tmp)
    if iscell(app.data(tmp(i)).name)
        name_str{i} = app.data(tmp(i)).name{1,1};
    else
        name_str{i} = app.data(tmp(i)).name;
    end
end
rem_idx = [];
if ~isnan(app.settings.channel_idx.ecg)
    rem_idx = [rem_idx ,find(contains(name_str,{'RR','Heartrate','heartrate'}))];
end

if ~isnan(app.settings.channel_idx.bldp)
    rem_idx = [rem_idx ,find(contains(name_str,{'BP', 'bp', 'Bloodpressure', 'bloodpressure', 'blood pressure', 'Blood Pressure'}))];
end

if ~isnan(app.settings.channel_idx.resp)
    rem_idx = [rem_idx ,find(contains(name_str,{'RESP', 'resp', 'Resp'}))];
end

if ~isnan(app.settings.channel_idx.co2)
    rem_idx = [rem_idx ,find(contains(name_str,{'co2', 'CO2'}))];
end

if ~isnan(app.settings.channel_idx.res_flow)
    rem_idx = [rem_idx ,find(contains(name_str,{'min. vent.', 'min vent', 'Minute Ventilation', 'minute ventilation'}))];
end

rem_idx = unique(rem_idx);
additional_sig = tmp;
additional_sig(rem_idx) = [];

rem_idx = [];
for i = 1 : length(app.stat_group)
    if isempty(app.stat_group(i).test)
        rem_idx = [rem_idx, i];
    end
end

app.stat_group(rem_idx) = [];

%% create results

plot_data = struct('name',[],'simple_name',[],'integral',[],'amplitude',[],'duration',[],...
                   'n_bursts',[], 'burstrate',[], 'burstincidence',[], ...
                   'mean_hr',[], 'add_sig',{});

for i = 1 : length(app.burst_ints)
    plot_data(i).name = app.burst_ints(i).name;%plot_data(i).name = app.burst_ints(sorted_idx(i)).name;
    plot_data(i).simple_name = simple_name(plot_data(i).name);
    tmp_idx_burst = app.burst_res.use_burst(:,1) & app.burst_res.use_burst(:,2) & app.burst_res.use_burst(:,i+2);
    dur = sum(diff(app.burst_ints(i).borders,1,2));%dur = sum(diff(app.burst_ints(sorted_idx(i)).borders,1,2));
    plot_data(i).integral = app.burst_res.burst_int(tmp_idx_burst);
    plot_data(i).amplitude = app.burst_res.burst_amp(tmp_idx_burst);
    plot_data(i).duration = app.burst_res.burst_dur(tmp_idx_burst);
    plot_data(i).latency =app.burst_res.burst_latency(tmp_idx_burst);
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

for i = 1 : length(additional_sig)
    [tmp_data,tmp_ts,~, ~] = current_signal(app, additional_sig(i));
    tmp_data (:,2) = tmp_ts(1) : tmp_ts(1): tmp_ts(2);
   
    for j = 1 : length(app.burst_ints)
        tmp_borders = app.burst_ints(j).borders;
        plot_data(j).add_sig{i} = tmp_data(tmp_data(:,2)>= tmp_borders(1) & tmp_data(:,2)<= tmp_borders(2),1);
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
tmp_cell_9 = {[]};
for i = 1:length(plot_data)
    namecell{i} = plot_data(i).simple_name;
    tmp_cell_1{i} = plot_data(i).integral;
    tmp_cell_2{i} = plot_data(i).amplitude;
    tmp_cell_3{i} = plot_data(i).duration;
    tmp_cell_9{i} = plot_data(i).latency;
    
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

plot_boxplot(tmp_cell_9,namecell,path,[file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) '_burst_analysis_burst_latency'],indx, 'Burst Latency')
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
plot_cell{14,1} = 'median latency [s]';
plot_cell{15,1} = 'mean latency [s]';
plot_cell{16,1} = 'std latency [s]';

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

if ~isnan(app.settings.channel_idx.res_flow)
    minvent_idx = size(plot_cell,1)+1;
    plot_cell{minvent_idx,1} = 'median minute ventilation';
    plot_cell{minvent_idx+1,1} = 'mean minute ventilation';
    plot_cell{minvent_idx+2,1} = 'std minute ventilation';
else
    minvent_idx = nan;
end

if ~isempty(additional_sig)
    add_idx = size(plot_cell,1)+1;
    for i  = 1 : length(additional_sig)
        if iscell(app.data(additional_sig(i)).name)
            tmp_name = app.data(additional_sig(i)).name{1,1};
        else
            tmp_name = app.data(additional_sig(i)).name;
        end

        plot_cell{add_idx+(i-1)*3,1} = ['median ' tmp_name];
        plot_cell{add_idx+(i-1)*3+1,1} = ['mean ' tmp_name];
        plot_cell{add_idx+(i-1)*3+2,1} = ['std ' tmp_name];
    end
else
    add_idx = nan;
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
    plot_cell{14,i+1} = median(plot_data(i).latency);
    plot_cell{15,i+1} = mean(plot_data(i).latency);
    plot_cell{16,i+1} = std(plot_data(i).latency);
   
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

    if  ~isnan(add_idx) 
        
        for j  = 1 : length(additional_sig)
            plot_cell{add_idx+(j-1)*3,i+1} = median(plot_data(i).add_sig{j});
            plot_cell{add_idx+(j-1)*3+1,i+1} = mean(plot_data(i).add_sig{j});
            plot_cell{add_idx+(j-1)*3+2,i+1} = std(plot_data(i).add_sig{j});
        end

    end

end

writecell(plot_cell,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],'FileType','spreadsheet','Sheet', 'general','Range','A1','WriteMode','replacefile');

app.lbl_working.Text = [num2str(round((1) / (1+length(app.stat_group))*100)) '% done'];
drawnow
rem_idx = [];
for i = 1 : length(additional_sig)
    if contains(app.data(additional_sig(i)).name,{'BURST INTEGRAL','BURST DURATION','BURST amplitude','BURST INTERBURST INTERVAL'})
        rem_idx = [rem_idx, i];
    end
end
additional_sig(rem_idx) = [];

for i =1:length(plot_data)
    plot_data(i).add_sig(rem_idx) = [];  
end

if ~isempty(app.stat_group)
    filename = [path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ];
    for i = 1 : length(app.stat_group)
        res = {[]};
        end_idx = 1;
        switch app.stat_group(i).paired
            case true
                
                switch app.stat_group(i).test
                    case 'two sample t test'
                        test_string = 'two sample t test';
                    case 'Wilcoxon rank sum test'
                        test_string = 'Wilcoxon rank sum test';
                end
                data1_name = plot_data(app.stat_group(i).intervals(1)).simple_name;
                data2_name = plot_data(app.stat_group(i).intervals(2)).simple_name;
                res{end_idx,1} = ['statistics group ' num2str(i)];
                res{end_idx,2} = data1_name; 
                res{end_idx,3} = ' vs. ';
                res{end_idx,4} = data2_name; 
                res{end_idx+1,1} = test_string;
                res{end_idx+1,2} = 'null hypothesis';
                res{end_idx+1,3} = 'p-value';
                end_idx = end_idx+1;
                int_idx = app.stat_group(i).intervals;
                plot_fields = {'E3', 'K3', 'Q3', 'W3', 'AC3', 'AJ3', 'AP3', 'AV3', 'BB3', 'BH3',...
                               'BO3', 'BU3', 'CA3', 'CG3', 'CM3', 'CS3', 'CY3', 'DE3', 'DK3', 'DQ3',...
                               'DW3', 'EC3', 'EJ3', 'EP3', 'EV3', 'FB3', 'FH3', 'FO3', 'FU3', 'GA3',...
                               'GG3', 'GM3', 'GS3', 'GY3', 'HE3', 'HK3', 'HQ3', 'HW3', 'IC3', 'IJ3'};
                plot_idx = 1;

                [res, end_idx] = paired_tests(plot_data(app.stat_group(i).intervals(1)).integral,...
                                              plot_data(app.stat_group(i).intervals(2)).integral,...
                                              data1_name, data2_name, res, end_idx,...
                                              plot_fields{plot_idx}, 'Burst Integral', i, app, filename);
                plot_idx = plot_idx+1;

                [res, end_idx] = paired_tests(plot_data(app.stat_group(i).intervals(1)).amplitude,...
                                              plot_data(app.stat_group(i).intervals(2)).amplitude,...
                                              data1_name, data2_name, res, end_idx,...
                                              plot_fields{plot_idx}, 'Burst Amplitude', i, app, filename);
                plot_idx = plot_idx+1;

                [res, end_idx] = paired_tests(plot_data(app.stat_group(i).intervals(1)).duration,...
                                              plot_data(app.stat_group(i).intervals(2)).duration,...
                                              data1_name, data2_name, res, end_idx,...
                                              plot_fields{plot_idx}, 'Burst Duration', i, app, filename);
                plot_idx = plot_idx+1;

               
                if ~isempty(app.hb_res)
                    [res, end_idx] = paired_tests(tmp_cell_4{app.stat_group(i).intervals(1)},...
                                                  tmp_cell_4{app.stat_group(i).intervals(2)},...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'RR-intervals', i, app, filename);
                    plot_idx = plot_idx+1;


                    
                    [res, end_idx] = paired_tests(tmp_cell_5{app.stat_group(i).intervals(1)},...
                                                  tmp_cell_5{app.stat_group(i).intervals(2)},...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'averaged heart rate', i, app, filename);
                    plot_idx = plot_idx+1;


                end

                if ~isnan(bp_idx)
                    
                    tmp_sig1 = [];
                    tmp_ts = app.data(sys_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig1 = [tmp_sig1, [app.data(sys_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)';
                                               app.data(mea_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)';
                                               app.data(dia_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)']];
                    end

                    tmp_sig2 = [];
                    tmp_ts = app.data(sys_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig2 = [tmp_sig2, [app.data(sys_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)';
                                               app.data(mea_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)';
                                               app.data(dia_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)']];
                    end

                    [res, end_idx] = paired_tests(tmp_sig1(1,:),...
                                                  tmp_sig2(1,:),...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'systolic BP', i, app, filename);
                    plot_idx = plot_idx+1;

                    [res, end_idx] = paired_tests(tmp_sig1(2,:),...
                                                  tmp_sig2(2,:),...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'mean BP', i, app, filename);
                    plot_idx = plot_idx+1;

                    [res, end_idx] = paired_tests(tmp_sig1(3,:),...
                                                  tmp_sig2(3,:),...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'diastolic BP', i, app, filename);
                    plot_idx = plot_idx+1;
                end
            
                if ~isnan(resp_idx)
                    
                    tmp_sig1 = [];
                    tmp_ts = app.data(resp_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig1 = [tmp_sig1, app.data(resp_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)'];
                    end

                    tmp_sig2 = [];
                    tmp_ts = app.data(resp_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig2 = [tmp_sig2, app.data(resp_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)'];
                    end

                    [res, end_idx] = paired_tests(tmp_sig1(1,:),...
                                                  tmp_sig2(1,:),...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'respiration rate ', i, app, filename);
                    plot_idx = plot_idx+1;

                end
               
%                 if ~isnan(baro_idx)
% 
%                 end
            
                if ~isnan(etco2_idx)  

                    tmp_sig1 = [];
                    tmp_ts = app.data(co2_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig1 = [tmp_sig1, app.data(co2_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)'];
                    end

                    tmp_sig2 = [];
                    tmp_ts = app.data(co2_ch_idx).ts(1);
                    for j = 1:size(app.burst_ints(int_idx(1)).borders,1)
                        tmp_sig2 = [tmp_sig2, app.data(co2_ch_idx).data(app.burst_ints(int_idx(1)).borders(j,1)/tmp_ts:app.burst_ints(int_idx(1)).borders(j,2)/tmp_ts)'];
                    end

                    [res, end_idx] = paired_tests(tmp_sig1(1,:),...
                                                  tmp_sig2(1,:),...
                                                  data1_name, data2_name, res, end_idx,...
                                                  plot_fields{plot_idx}, 'etCO2', i, app, filename);
                    plot_idx = plot_idx+1;

                end
            
                if  ~isnan(add_idx) 

                    use_sigs = [];
                    for j = 1:length(additional_sig)
                        if ismember(additional_sig(j), app.stat_group(i).signals)
                            use_sigs = [use_sigs, j];
                        end
                    end

                    for j = 1:length(use_sigs)
                        if iscell(app.data(additional_sig(j)).name)
                            signame = app.data(additional_sig(j)).name{1,1};
                        else
                            signame = app.data(additional_sig(j)).name;
                        end
                        
                        [res, end_idx] = paired_tests(plot_data(int_idx(1)).add_sig{use_sigs(j)},...
                                                      plot_data(int_idx(2)).add_sig{use_sigs(j)},...
                                                      data1_name, data2_name, res, end_idx,...
                                                      plot_fields{plot_idx}, signame, i, app, filename);
                        plot_idx = plot_idx+1;

                    end

                end

                writecell(res, filename, 'FileType','spreadsheet','Sheet', ['statistics group ' num2str(i) ],'Range','A1')
    %% anova
            case false
                names = [];
                namecell= {[]};
                int_idx = app.stat_group(i).intervals;
                for k = 1 : length(int_idx)
                    if k == 1
                        names = [simple_name(app.burst_ints(int_idx(k)).name)];
                    else
                        names = [names ' vs. ' simple_name(app.burst_ints(int_idx(k)).name)];
                    end
                    namecell{k} = simple_name(app.burst_ints(int_idx(k)).name);
                end       

                if strcmp(app.stat_group(i).posthoc , 'none')
                    res{1,1} = ['statistics group ' num2str(i)   ];
                    res{1,2} = names; 
                    res{2,1} = app.stat_group(i).test;
                    res{2,2} = 'p-value';
                    
                else
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

                plot_fields = {'E1','H1'; 'E21', 'H23'; 'E41','H45'; 'E61','H67'; 'E81', 'H89'; ...
                               'E101', 'H111'; 'E121', 'H133'; 'E141', 'H155'; 'E161', 'H177'; 'E181', 'H199';...
                               'E201', 'H221'; 'E221', 'H243'; 'E241', 'H265'; 'E261', 'H287'; 'E281', 'H309';...
                               'E301', 'H331'; 'E321', 'H353'; 'E341', 'H375'; 'E361', 'H397'; 'E381', 'H419';...
                               'E401', 'H441'; 'E421', 'H463'; 'E441', 'H485'; 'E461', 'H507'; 'E481', 'H529';...
                               'E501', 'H551'; 'E521', 'H573'; 'E541', 'H595'; 'E561', 'H617'; 'E581', 'H639';...
                               'E601', 'H661'; 'E621', 'H683'; 'E641', 'H663'; 'E661', 'H685'; 'E681', 'H707'};
                num_plot = 0;
                end_idx = 2;

                vals = [];
                group = string([]);
                tmp_cell = {[]};
                for k = 1 : length(int_idx)
                    vals = [vals,plot_data(int_idx(k)).integral];
                    group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                    tmp_cell{k} = plot_data(int_idx(k)).integral;
                end                
                [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'Burst integral', i, plot_fields, end_idx, num_plot, app, filename, namecell);
                
                vals = [];
                group = string([]);
                tmp_cell = {[]};
                for k = 1 : length(int_idx)
                    vals = [vals,plot_data(int_idx(k)).amplitude];
                    group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                    tmp_cell{k} = plot_data(int_idx(k)).amplitude;
                end                
                [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'Burst Amplitude', i, plot_fields, end_idx, num_plot, app, filename, namecell);
                
                vals = [];
                group = string([]);
                tmp_cell = {[]};
                for k = 1 : length(int_idx)
                    vals = [vals,plot_data(int_idx(k)).duration'];
                    group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                    tmp_cell{k} = plot_data(int_idx(k)).duration;
                end                
                [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'Burst Duration', i, plot_fields, end_idx, num_plot, app, filename, namecell);
               
                if ~isempty(app.hb_res)
                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_cell_4{int_idx(k)}'];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_cell_4{int_idx(k)};
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'RR-intervals', i, plot_fields, end_idx, num_plot, app, filename, namecell);
    
                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_cell_5{int_idx(k)}'];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_cell_5{int_idx(k)};
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'averaged heart rate', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                end

                if ~isnan(bp_idx)
                    tmp_data = {};

                    for j = 1: length(int_idx)

                        tmp_sig = [];
                        tmp_ts = app.data(sys_ch_idx).ts(1);
                        for k = 1:size(app.burst_ints(int_idx(j)).borders,1)
                            tmp_sig = [tmp_sig, [app.data(sys_ch_idx).data(app.burst_ints(int_idx(1)).borders(k,1)/tmp_ts:app.burst_ints(int_idx(j)).borders(k,2)/tmp_ts)';
                                                   app.data(mea_ch_idx).data(app.burst_ints(int_idx(1)).borders(k,1)/tmp_ts:app.burst_ints(int_idx(j)).borders(k,2)/tmp_ts)';
                                                   app.data(dia_ch_idx).data(app.burst_ints(int_idx(1)).borders(k,1)/tmp_ts:app.burst_ints(int_idx(j)).borders(k,2)/tmp_ts)']];
                        end
                        tmp_data{j} = tmp_sig;
                    end
                    
                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_data{k}(1,:)];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_data{k}(1,:);
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'systolic BP', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_data{k}(2,:)];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_data{k}(2,:);
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'mean BP', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_data{k}(3,:)];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_data{k}(3,:);
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'diastolic BP', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                end
            
                if ~isnan(resp_idx)

                    tmp_data = {};

                    for j = 1: length(int_idx)

                        tmp_sig = [];
                        tmp_ts = app.data(resp_ch_idx).ts(1);
                        for k = 1:size(app.burst_ints(int_idx(j)).borders,1)
                            tmp_sig = [tmp_sig, app.data(resp_ch_idx).data(app.burst_ints(int_idx(1)).borders(k,1)/tmp_ts:app.burst_ints(int_idx(j)).borders(k,2)/tmp_ts)];
                        end
                        tmp_data{j} = tmp_sig;
                    end

                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_data{k}'];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_data{k};
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'respiration rate', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                end
               
%                 if ~isnan(baro_idx)
% 
%                 end
            
                if ~isnan(etco2_idx)  
                    
                    tmp_data = {};

                    for j = 1: length(int_idx)

                        tmp_sig = [];
                        tmp_ts = app.data(co2_ch_idx).ts(1);
                        for k = 1:size(app.burst_ints(int_idx(j)).borders,1)
                            tmp_sig = [tmp_sig, app.data(co2_ch_idx).data(app.burst_ints(int_idx(1)).borders(k,1)/tmp_ts:app.burst_ints(int_idx(j)).borders(k,2)/tmp_ts)'];
                        end
                        tmp_data{j} = tmp_sig;
                    end

                    vals = [];
                    group = string([]);
                    tmp_cell = {[]};
                    for k = 1 : length(int_idx)
                        vals = [vals,tmp_data{k}];
                        group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                        tmp_cell{k} = tmp_data{k};
                    end                
                    [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, 'etCO2', i, plot_fields, end_idx, num_plot, app, filename, namecell);

                end
            
                if  ~isnan(add_idx) 

                     use_sigs = [];
                    for j = 1:length(additional_sig)
                        if ismember(additional_sig(j), app.stat_group(i).signals)
                            use_sigs = [use_sigs, j];
                        end
                    end

                    for j = 1:length(use_sigs)
                        if iscell(app.data(additional_sig(j)).name)
                            signame = app.data(additional_sig(j)).name{1,1};
                        else
                            signame = app.data(additional_sig(j)).name;
                        end
                        
                        vals = [];
                        group = string([]);
                        tmp_cell = {[]};
                        for k = 1 : length(int_idx)
                            vals = [vals,plot_data(int_idx(k)).add_sig{use_sigs(j)}'];
                            group(end+1:length(vals)) = string(simple_name(app.burst_ints(int_idx(k)).name));
                            tmp_cell{k} = plot_data(int_idx(k)).add_sig{use_sigs(j)};
                        end                
                        [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, signame, i, plot_fields, end_idx, num_plot, app, filename, namecell);

                    end


                end

                writecell(res,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_burst_results.xls' ],'FileType','spreadsheet','Sheet', ['statistics group ' num2str(i) ],'Range','A1')
            
        end
        app.lbl_working.Text = [num2str(round((1+i) / (1+length(app.stat_group))*100)) '% done'];
        drawnow
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
            h.Visible = 'on';
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


function [res, end_idx] = paired_tests(data1, data2, data1_name, data2_name, res,end_idx, plot_field, title, stat_group,app, filename)
print_idx = end_idx+1; 
switch app.stat_group(stat_group).test
    case 'two sample t test'
        [h,p] =  ttest2(data1, data2);
    case 'Wilcoxon rank sum test'
        [p,h] = ranksum(data1, data2);
end


res{print_idx,1} = title;
if islogical(h)
    h = int8(h);
end
res{print_idx,2} = h;
res{print_idx,3} = p;
stat_plot_cell = {data1, data2};

plot_boxplot_xls(stat_plot_cell, {data1_name, data2_name}, ...
                 title, filename, ['statistics group ' num2str(stat_group) ], h, p, plot_field)
end_idx = print_idx;
end


function [res, end_idx, num_plot] = anova_tests(vals, res, group, tmp_cell, title, stat_group, plot_fields, end_idx, num_plot, app, filename, namecell)
num_plot = num_plot +1;
end_idx = end_idx +1;

switch app.stat_group(stat_group).test
    case 'anova'
        [p,~,stats] = anova1(vals,group,"off");
    case 'anocova'
        %[h,atab,ctab,stats] = aoctool(x,y,group,alpha,xname,yname,gname,"off")
end
                
if strcmp(app.stat_group(stat_group).posthoc , 'none')
    
    res{end_idx,1} = title;
    res{end_idx,2} = p;
    plot_boxplot_xls(tmp_cell,namecell , title, filename, ['statistics group ' num2str(stat_group) ], false, p, plot_fields{num_plot,1})

else
    c = multcompare(stats, 'CriticalValueType', app.stat_group(stat_group).posthoc);
    p = c(c(:,6)<=0.05,[1,2,6]);
    
    res{end_idx,1} = title;
    end_idx = end_idx+1;
    cc = num2cell(c);
    for k = 1:size(c,1)
        cc{k,1} = namecell{c(k,1)};
        cc{k,2} = namecell{c(k,2)};
    end
    res(end_idx:end_idx+size(c,1)-1,1:6) = cc; 
    end_idx = end_idx+size(c,1)-1;
    plot_boxplot_phc_xls(tmp_cell,namecell , title, filename, ['statistics group ' num2str(stat_group) ],  p, plot_fields{num_plot,2})
    
end

end