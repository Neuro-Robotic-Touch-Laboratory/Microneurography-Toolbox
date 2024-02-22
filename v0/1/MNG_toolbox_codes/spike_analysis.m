function spike_res = spike_analysis(app)
%SPIKE_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here


MIN_CLUSTER=10;


calc_flag = true;

if ~isempty(app.spike_res)
    spike_res = app.spike_res;
    if ~spike_res.analysis.sorting && ~app.chkbx_spike_sorting.Value
        calc_flag = false;
    end
    if (spike_res.analysis.sorting && app.chkbx_spike_sorting.Value) && isequal(spike_res.analysis.dips,  [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value])
        calc_flag = false;
    end
    
else
    spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
    spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
end

if calc_flag
    data = app.data(app.settings.channel_idx.msna).data(app.settings.interval(1,1)/app.data(app.settings.channel_idx.msna).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.msna).ts(1));

    
    data = data';
    writable_folder = GetWritableFolder;
    save ([writable_folder '\temp.mat'], "data") %save ('temp.mat', "data")
    
    SimpleSpikesSorting = app.chkbx_spike_sorting.Value;
    
    Get_spikes_folder([writable_folder '\temp.mat'])% Get_spikes([writable_folder '\temp.mat']) % add sample rate
    spikes=load([writable_folder '\temp_spikes.mat'],'spikes','index','threshold');
    

    if app.chkbx_spike_sorting.Value
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
    spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
    spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
    delete([writable_folder '\temp_spikes.mat'])
    delete([writable_folder '\temp.mat'])
    spike_res = sel_spikes(spike_res, app.data(app.settings.channel_idx.msna).ts, app.burst_ints); 
end

%% spike analysis
if ~spike_res.analysis.spike
   
    for k = 1:length(app.popup_int_spike.Items)

        radX = {[]};
        SpikeTimeX = {[]};
        first = {[]};
        last = {[]};
        totXrad_mean = {[]};
        totXrad_min = {[]};
        totXrad_max = {[]};
        bp_fr_base = {[]};
        firstx10 = {[]};
        avgx10 = {[]};
        lastx10 = {[]};
        shape = {[]};
        n_spikes = nan (max(spike_res.cluster),1);
        tmp_rr = [nan(10,1);diff(app.hb_res.t_events(app.hb_res.use_beats(:,k))); nan(10,1)];
        calcmat = zeros(length(tmp_rr),20);
        for j =1:9
            calcmat(:,1:10-j) = calcmat(:,1:10-j)+circshift(tmp_rr,j);
        end
        for j = 1:10
            calcmat(:,10+j:20) = calcmat(:,10+j:20)+circshift(tmp_rr,-j);
        end
        if k == 1
            dur = diff(app.settings.interval(1,:));
        else
            dur = sum (diff(app.burst_ints(k-1).borders,1,2));
        end

        for i = 1: max(spike_res.cluster)
            if k == 1
                tmp_spk_idx = spike_res.cluster == i;
            else
                tmp_spk_idx = (spike_res.cluster == i) & spike_res.use_spikes(:,k)'; 
            end
            n_spikes(i) = sum(tmp_spk_idx);
            FR_tot = sum(tmp_spk_idx)/dur;
            shape{i} = mean(spike_res.spikes(:,tmp_spk_idx),2);
            [radX{i}, SpikeTimeX{i}, first{i}, last{i}] = computeRadFiring(app.hb_res.t_events(app.hb_res.use_beats(:,k)),spike_res.spike_ts(:, tmp_spk_idx)/1000);
           
            for j = 1 : length(first{i})
                tmp_idx = find(spike_res.spike_ts/1000 >= first{i}(j,1),1);
                if ~isempty(tmp_idx)
                    first{i}(j,2) = tmp_idx;
                else
                    first{i}(j,2) = nan;
                end
                tmp_idx = find(spike_res.spike_ts/1000 >=  last{i}(j,1),1);
                if ~isempty(tmp_idx)
                    last{i}(j,2) = tmp_idx;
                else
                    last{i}(j,2) = nan;
                end
            end
            dt_cycle = nan(3,length (SpikeTimeX{i}));
            for j = 1 : length (SpikeTimeX{i})
                if ~isempty(SpikeTimeX{i}{1,j})
                    dt_cycle(1,j) = SpikeTimeX{i}{1,j}(1);
                    dt_cycle(2,j) = mean(SpikeTimeX{i}{1,j});
                    dt_cycle(3,j) = SpikeTimeX{i}{1,j}(end);
                end
            end
            dt_cycle = [nan(3,10), dt_cycle, nan(3,10)];

            firstx10{i} = calcmat(:,1:10)+dt_cycle(1,:)';
            firstx10{i}(:,11:20) = calcmat(:,11:20)-dt_cycle(1,:)';
            firstx10{i}(end-9:end,:) = [];
            firstx10{i}(1:10,:) = [];
            avgx10{i} = calcmat(:,1:10)+dt_cycle(2,:)';
            avgx10{i}(:,11:20) = calcmat(:,11:20)-dt_cycle(2,:)';
            avgx10{i}(end-9:end,:) = [];
            avgx10{i}(1:10,:) = [];
            lastx10{i} = calcmat(:,1:10)+dt_cycle(3,:)';
            lastx10{i}(:,11:20) = calcmat(:,11:20)-dt_cycle(3,:)';
            lastx10{i}(end-9:end,:) = [];
            lastx10{i}(1:10,:) = [];


            [totXrad_mean{i}, totXrad_min{i}, totXrad_max{i}]= ComputeLatencyCell(radX{i});
            dt_cycle = diff(app.hb_res.t_events(app.hb_res.use_beats(:,k)));
            
            %FR_tot=length(t_SPIKES)/t(end);
            tmp_spk_num = nan(1,length(radX{i}));
            for l=1:length(radX{i})
                tmp_spk_num(l) = length(radX{i}{1,l});

            end      
            
            
            r = circ_r(totXrad_min{i});
            phi = circ_mean(totXrad_min{i});
            if phi<0
                phi=2*pi-abs(phi);
            end
            circFIRST=mean(phi);
            plvFIRST=mean(r);            

            r = circ_r(totXrad_mean{i});
            phi = circ_mean(totXrad_mean{i});
            if phi<0
                phi = 2*pi-abs(phi);
            end
            circAVG=mean(phi);
            plvAVG=mean(r);

            r = circ_r(totXrad_max{i});
            phi = circ_mean(totXrad_max{i});
            if phi<0
                phi=2*pi-abs(phi);
            end
            
            circLAST=mean(phi);
            plvLAST=mean(r);
            phase2first{i} = [FR_tot, mean(totXrad_mean{i}), mean(totXrad_min{i}), mean(totXrad_max{i}), circAVG, circFIRST, circLAST, plvAVG, plvFIRST, plvLAST];

        end


        spike_res.spike(k).radX = radX;
        spike_res.spike(k).SpikeTimeX = SpikeTimeX; 
        spike_res.spike(k).first = first;
        spike_res.spike(k).last = last;
        spike_res.spike(k).totXrad_mean = totXrad_mean;
        spike_res.spike(k).totXrad_min = totXrad_min;
        spike_res.spike(k).totXrad_max = totXrad_max;
        spike_res.spike(k).analysis.spike = true; 
%         spike_res.spike(k).bp_fr_base = bp_fr_base;
%         spike_res.spike(k).firstx10_all = firstx10_all;
%         spike_res.spike(k).avgx10_all = avgx10_all;
%         spike_res.spike(k).lastx10_all = lastx10_all;
        spike_res.spike(k).firstx10 = firstx10;
        spike_res.spike(k).avgx10 = avgx10;
        spike_res.spike(k).lastx10 = lastx10;
        spike_res.spike(k).phase2first = phase2first;
        spike_res.spike(k).n_spikes = n_spikes;
        spike_res.spike(k).shape = shape;
    end
end

end










