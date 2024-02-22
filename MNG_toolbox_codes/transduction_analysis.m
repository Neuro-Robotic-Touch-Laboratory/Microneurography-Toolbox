function spike_res = transduction_analysis (app)


MIN_CLUSTER=10;
% MSNA = app.data(app.settings.channel_idx.msna).data;
% data = app.data(app.settings.channel_idx.msna).data(app.settings.interval(1,1)/app.data(app.settings.channel_idx.msna).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.msna).ts(1));

calc_flag = true;

if ~isempty(app.spike_res)
    spike_res = app.spike_res;

    if (spike_res.analysis.sorting == app.chkbx_sorting_trans.Value) ...
        && isequal(spike_res.analysis.dips,  [app.edt_trans_dip1.Value, app.edt_trans_peak.Value, app.edt_trans_dip2.Value]) ...
        && (spike_res.analysis.waveclus == app.chkbx_trans_clustering.Value...
        && (spike_res.analysis.inverse == app.chbx_invert_trans.Value))
        calc_flag = false;
    end
    
    if calc_flag
         
        spike_res.analysis.sorting = app.chkbx_sorting_trans.Value;
        spike_res.analysis.waveclus = app.chkbx_trans_clustering.Value;
        spike_res.analysis.dips = [app.edt_trans_dip1.Value, app.edt_trans_peak.Value, app.edt_trans_dip2.Value];
        spike_res.analysis.inverse = app.chbx_invert_trans.Value;
        spike_res.analysis.spike = false;
        spike_res.analysis.trans = false;
        spike_res.spikes = [];
        spike_res.spike_idx = [];
        spike_res.spike_ts = [];
        spike_res.cluster = [];
        spike_res.extremes = [];
        spike_res.use_spikes = [];
        spike_res.spike = [];
        spike_res.rate_ylim = [];
    end
else
    spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
    spike_res.analysis.waveclus = app.chkbx_spike_clustering.Value;
    spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
    spike_res.analysis.inverse = app.chbx_invert_spike.Value;
end

if calc_flag
    %     data = app.data(app.settings.channel_idx.msna).data(app.settings.interval(1,1)/app.data(app.settings.channel_idx.msna).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.msna).ts(1));
    [data,ts,~,~] = current_signal(app, app.settings.channel_idx.msna);
    if app.chbx_invert_spike.Value
        data = data * (-1);
    end
    
    data = data';
    writable_folder = GetWritableFolder;
    sr = 1/ts(1);
    save ([writable_folder '\temp.mat'], "data","sr") %save ('temp.mat', "data")
    
    SimpleSpikesSorting = app.chkbx_spike_sorting.Value;
    
    Get_spikes_folder([writable_folder '\temp.mat'])% Get_spikes([writable_folder '\temp.mat']) % add sample rate
    if app.chkbx_trans_clustering.Value
        app.edt_trans_file.Value = [writable_folder '\temp.mat'];
        clipboard('copy',[writable_folder '\temp.mat'])
        h = wave_clus;
        waitfor(h)
    end  

    if app.chkbx_trans_clustering.Value 
        cluster=load([writable_folder '\times_temp.mat'],'cluster_class');
        [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
                =comBIN_wave_clus(cluster,data,sr);
    else
        spikes=load([writable_folder '\temp_spikes.mat'],'spikes','index','threshold');
        if app.chkbx_sorting_trans.Value
            %%speed up with vector operation
            [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
                =comBIN12_v03(spikes,[writable_folder '\temp.mat'],app.edt_dip1.Value,app.edt_peak.Value,app.edt_dip2.Value);
    %         [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
    %             =comBIN12_v03(spikes,'temp',app.edt_dip1.Value,app.edt_peak.Value,app.edt_dip2.Value);
            %%
        else
            [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
                =comBIN12_v03(spikes,[writable_folder '\temp.mat'],1,1,1);
    %         [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
    %             =comBIN12_v03(spikes,'temp',1,1,1);
        end
    end
    spike_res.analysis.sorting = app.chkbx_sorting_trans.Value;
    spike_res.analysis.dips = [app.edt_trans_dip1.Value, app.edt_trans_peak.Value, app.edt_trans_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
    spike_res = sel_spikes(spike_res, app.data(app.settings.channel_idx.msna).ts, app.burst_ints); 
    list = dir(writable_folder);
    list(1:2) = [];
    for i = 1: length(list)
        delete([list(i).folder '\' list(i).name])
    end
end
%%

if ~spike_res.analysis.trans
    
    int_names = app.popup_int_transduction.Items;
    ints_idxs =nan(size(int_names));
    all_ints = {[]};
    for i = 1: length(app.burst_ints)
        all_ints{i} = app.burst_ints(i).name;
    end
    borders = nan(length(int_names),2);
    for i = 2: length(int_names)
        ints_idxs(i) = find (contains(all_ints, int_names{i}));
        borders(i,:) = app.burst_ints(ints_idxs(i)).borders;
    end
    %%%%%%
%     basic_ints =find(vertcat(app.burst_ints.type) == 1);                %get intervals
%     borders = nan(length(basic_ints)+1,2);                              %from popup items 
%     int_names = {'full interval'};                                      %
%        
%     for i = 1: length(basic_ints)                                       %
%         borders(i+1,:) = app.burst_ints(basic_ints(i)).borders;         %
%         int_names{i+1} =app.burst_ints(basic_ints(i)).name;             %
%     end

    [data_msna,ts_msna,~, ~] = current_signal(app, app.settings.channel_idx.msna);
    if ~isnan(app.settings.channel_idx.bldp)
        [data_bp,ts_bp,~, ~] = current_signal(app, app.settings.channel_idx.bldp);
    else
        data_bp = nan;
        ts_bp = nan;
    end
    [~,ts_ecg,~, ~] = current_signal(app, app.settings.channel_idx.ecg);
    xxxR=-10:1:10;
    xxxR=nonzeros(xxxR);    
    clus = unique(spike_res.cluster);
    for i = 1 : length(int_names)

        for j = 1 : length(clus)
            if isnan(ints_idxs(i))
                tmp_idx = spike_res.cluster ==clus(j);
                t_spikes = spike_res.spike_ts (tmp_idx);
                data = data_msna';
                t_events_ecg = app.hb_res.t_events;
%                 if ~isnan(data_bp)
                    bpValues = data_bp;
                    foot_idx = app.bp_res.foot_idx;
                    t_bpFoot = foot_idx*ts_bp(1);
%                 else
%                     bpValues = nan;
%                     foot_idx = nan;
%                     t_bpFoot = nan;
%                 end

            else
                tmp_idx = spike_res.use_spikes(:,ints_idxs(i)+1) & (spike_res.cluster == clus(j))';%%% change for long intervals
                t_spikes = spike_res.spike_ts(tmp_idx) - (borders(i,1)-ts_msna(1))*1000;
                data = data_msna(borders(i,1)/ts_msna(1):borders(i,2)/ts_msna(1));
                t_events_ecg = app.hb_res.t_events(app.hb_res.use_beats(:,i))-borders(i,1) +ts_ecg(1);
%                 if ~isnan(data_bp)
                    bpValues = data_bp(borders(i,1)/ts_bp(1):borders(i,2)/ts_bp(1));
                    foot_idx = app.bp_res.foot_idx((app.bp_res.foot_idx >= borders (i,1)/ts_bp(1)) & (app.bp_res.foot_idx <= borders (i,2)/ts_bp(1))) -round(borders (i,1)/ts_bp(1));
                    t_bpFoot = foot_idx*ts_bp(1);
%                 else
%                     bpValues = nan;
%                     foot_idx = nan;
%                     t_bpFoot = nan;
%                 end
            end
            
            spike_res.transduction(i).cluster(j).mean_shape = mean(spike_res.spikes(:,tmp_idx),2);
            spike_res.transduction(i).cluster(j).n_spikes = sum(tmp_idx);
            
            [RRR7, RRR8 ] = transduction_calc(data, t_events_ecg, t_bpFoot, bpValues', t_spikes, foot_idx); %%%% adapt from GDAplot_DynamicIndexes_v26
            spike_res.transduction(i).cluster(j).phase_bp = [xxxR(11:20),RRR7(11:20,3)];
            spike_res.transduction(i).cluster(j).phase_rr = [xxxR(11:20),RRR7(11:20,4)];
            spike_res.transduction(i).cluster(j).latency_bp = [xxxR,RRR7(:,1)];
            spike_res.transduction(i).cluster(j).latency_rr = [xxxR,RRR7(:,2)];
            spike_res.transduction(i).cluster(j).fr_bp = [xxxR,RRR8(:,1)];
            spike_res.transduction(i).cluster(j).fr_rr = [xxxR,RRR8(:,2)];
        end
    end
    spike_res.analysis.trans = true;
end

end
