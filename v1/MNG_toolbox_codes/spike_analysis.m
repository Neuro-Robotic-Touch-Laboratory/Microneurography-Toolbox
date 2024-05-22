function spike_res = spike_analysis(app, update)
%SPIKE_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    update = false;
end
MIN_CLUSTER=10;


calc_flag = true;

% if ~isempty(app.spike_res)
%     spike_res = app.spike_res;
%     
%     if (spike_res.analysis.sorting == app.chkbx_spike_sorting.Value) ...
%         && isequal(spike_res.analysis.dips,  [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value]) ...
%         && (spike_res.analysis.waveclus == app.chkbx_spike_clustering.Value)
%         calc_flag = false;
%     end
%     
%     if calc_flag
%          
%         spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
%         spike_res.analysis.waveclus = app.chkbx_spike_clustering.Value;
%         spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
% %         spike_res.analysis.inverse = app.chbx_invert_spike.Value;
%         spike_res.analysis.spike = false;
%         spike_res.analysis.trans = false;
%         spike_res.spikes = [];
%         spike_res.spike_idx  = [];
%         spike_res.spike_ts  = [];
%         spike_res.cluster  = [];
%         spike_res.extremes  = [];
%         spike_res.use_spikes  = [];
%         spike_res.spike  = [];
%         spike_res.rate_ylim  = [];
%     end
% 
% else
    spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
    spike_res.analysis.waveclus = app.chkbx_spike_clustering.Value;
    spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
    spike_res.spikes = [];
    spike_res.spike_idx  = [];      %%
    spike_res.spike_ts  = [];       %%
    spike_res.cluster  = [];        %%
    spike_res.extremes  = [];       %%
    spike_res.use_spikes  = [];     %%
    spike_res.spike  = [];          %%
    spike_res.rate_ylim  = [];      %%
%     spike_res.analysis.inverse = app.chbx_invert_spike.Value;

% end

if calc_flag
    [data,ts,~,~] = current_signal(app, app.settings.channel_idx.msna);
    
    data = data'; 

    sr = 1/ts(1); % rem

    SimpleSpikesSorting = app.chkbx_spike_sorting.Value;
    
    if update
        threshold = app.spike_res.det_res.threshold;
        spikes = app.spike_res.det_res.spikes;
        spk_pos = app.spike_res.det_res.spk_pos;
        index = app.spike_res.det_res.index;
        use_beats = app.spike_res.use_spikes(:,1);
        par = app.settings.par;
    else
        [threshold, index, par, spikes, spk_pos] = detect_spikes(data, app.settings.par.sr, app.settings.par);
        app.settings.par = par;
    end

    disp([num2str(length(index)) 'length index']) %
    disp([num2str(size(spikes)) 'size spikes']) %
    if app.chkbx_spike_clustering.Value
        disp([num2str(length(index)) 'length index2']) %
        disp([num2str(size(spikes)) 'size spikes2']) %
        app.wc_app = cluster_app(app,index, spikes, par, threshold);
        waitfor(app.wc_app)
        cluster.cluster_class = app.settings.tempclus.cluster_class;
        app.settings.tempclus = [];
    end    
    
   
    if app.chkbx_spike_clustering.Value 
%         cluster=load([writable_folder '\times_temp.mat'],'cluster_class');
        [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
                =comBIN_wave_clus(cluster,data,sr);
    else
%         spikes=load([writable_folder '\temp_spikes.mat'],'spikes','index','threshold'); % rem
        spks =  struct('spikes',spikes,'index',index,'threshold',threshold);
        
        if app.chkbx_spike_sorting.Value
            %%speed up with vector operation
%             [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ... % rem
%                 =comBIN12_v03(spikes,[writable_folder '\temp.mat'],app.edt_dip1.Value,app.edt_peak.Value,app.edt_dip2.Value);                                           % rem
            [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts, spike_res.cluster, spike_res.extremes(:,2),  spike_res.extremes(:,1), spike_res.extremes(:,3) ]...
            =comBIN(spks,data,sr, spk_pos,app.edt_dip1.Value,app.edt_peak.Value,app.edt_dip2.Value);
            %%
        else
%             [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts,spike_res.cluster,spike_res.extremes(:,2), spike_res.extremes(:,1), spike_res.extremes(:,3)] ...
%                 =comBIN12_v03(spikes,[writable_folder '\temp.mat'],1,1,1);
            [spike_res.spikes, spike_res.spike_idx, spike_res.spike_ts, spike_res.cluster, spike_res.extremes(:,2),  spike_res.extremes(:,1), spike_res.extremes(:,3) ]...
            =comBIN(spks,data,sr,spk_pos,1,1,1);
        end
    end
    spike_res.analysis.sorting = app.chkbx_spike_sorting.Value;
    spike_res.analysis.dips = [app.edt_dip1.Value, app.edt_peak.Value, app.edt_dip2.Value];
    spike_res.analysis.spike = false;
    spike_res.analysis.trans = false;
    spike_res = sel_spikes(spike_res, app.data(app.settings.channel_idx.msna).ts, app.burst_ints);
    if update 
        spike_res.use_spikes(:,1) = app.spike_res.use_spikes(:,1);
    end
    spike_res.det_res.threshold = threshold;
    spike_res.det_res.spk_pos = spk_pos;
    spike_res.det_res.spikes = spikes;
    spike_res.det_res.index = index;
%     list = dir(writable_folder);
%     list(1:2) = [];
%     for i = 1: length(list)
%         delete([list(i).folder '\' list(i).name])
%     end
end

%% spike analysis
if ~spike_res.analysis.spike 
    spike_res.spike = [];
    for k = 1:length(app.popup_int_spike.Items) %% check  basic intor also derived ints

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
        m_fr = [];
        n_spikes = nan (max(spike_res.cluster),1);
        if ~isnan(app.settings.channel_idx.ecg)
            usebeats = app.hb_res.use_beats(:,k);
            if sum(usebeats) <=1
                tmp_idx = find(usebeats);
                if tmp_idx >1
                    usebeats(tmp_idx-1) = true;
                end
                if tmp_idx < length(usebeats)
                    usebeats(tmp_idx+1) = true;
                end

            end

            tmp_rr = [nan(10,1);diff(app.hb_res.t_events(usebeats)); nan(10,1)];
            % tmp_rr = [nan(10,1);diff(app.hb_res.t_events(app.hb_res.use_beats(:,k))); nan(10,1)];
            calcmat = zeros(length(tmp_rr),20);
            for j =1:9
                calcmat(:,1:10-j) = calcmat(:,1:10-j)+circshift(tmp_rr,j);
            end
            for j = 1:10
                calcmat(:,10+j:20) = calcmat(:,10+j:20)+circshift(tmp_rr,-j);
            end
        end
        if k == 1
            dur = diff(app.settings.interval(1,:));
        else
            dur = sum (diff(app.burst_ints(k-1).borders,1,2)); %% works with derived ??
        end

        for i = 1: max(spike_res.cluster)
            if k == 1
                tmp_spk_idx = spike_res.cluster == i;
            else
                tmp_spk_idx = (spike_res.cluster == i) & spike_res.use_spikes(:,k)'; 
            end
            n_spikes(i) = sum(tmp_spk_idx);
            FR_tot = sum(tmp_spk_idx)/dur;
            m_fr(i) = FR_tot;
            shape{i} = mean(spike_res.spikes(:,tmp_spk_idx),2);
            if ~isnan(app.settings.channel_idx.ecg)
                %disp(num2str([k,i]))
                [radX{i}, SpikeTimeX{i}, first{i}, last{i}] = computeRadFiring(app.hb_res.t_events(usebeats),spike_res.spike_ts(:, tmp_spk_idx)/1000); 
                % [radX{i}, SpikeTimeX{i}, first{i}, last{i}] = computeRadFiring(app.hb_res.t_events(app.hb_res.use_beats(:,k)),spike_res.spike_ts(:, tmp_spk_idx)/1000);
               
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
                dt_cycle = diff(app.hb_res.t_events(usebeats));
                % dt_cycle = diff(app.hb_res.t_events(app.hb_res.use_beats(:,k)));
            
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
        end
        if ~isnan(app.settings.channel_idx.ecg)
            spike_res.spike(k).radX = radX;
            spike_res.spike(k).SpikeTimeX = SpikeTimeX; 
            spike_res.spike(k).first = first;
            spike_res.spike(k).last = last;
            spike_res.spike(k).totXrad_mean = totXrad_mean;
            spike_res.spike(k).totXrad_min = totXrad_min;
            spike_res.spike(k).totXrad_max = totXrad_max;
            spike_res.spike(k).firstx10 = firstx10;
            spike_res.spike(k).avgx10 = avgx10;
            spike_res.spike(k).lastx10 = lastx10;
            spike_res.spike(k).phase2first = phase2first;
            
        end
        spike_res.spike(k).analysis.spike = true; 
%         spike_res.spike(k).bp_fr_base = bp_fr_base;
%         spike_res.spike(k).firstx10_all = firstx10_all;
%         spike_res.spike(k).avgx10_all = avgx10_all;
%         spike_res.spike(k).lastx10_all = lastx10_all;
        spike_res.spike(k).m_fr = m_fr;
        spike_res.spike(k).n_spikes = n_spikes;
        spike_res.spike(k).shape = shape;
    end
end

end










