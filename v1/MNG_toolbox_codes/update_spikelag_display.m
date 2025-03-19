function update_spikelag_display(app, status)

edges = -0.5 :0.005:3;
switch status
    case 'new'
        spks = app.spike_res.spike_lag.spks;
        num_beats = app.edt_num_beat.Value;
        mov_mean = app.edt_lag_movmean.Value;
        cluster = find(strcmp(app.popup_cluster_lag.Value,app.popup_cluster_lag.Items));         
        int_idx = find(strcmp(app.popup_int_spike.Value,app.popup_int_spike.Items));
        res_sm= nan(length(edges)-1,length(spks)-num_beats+2);
        peak = nan(1,length(spks)-num_beats+1);
        int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
        tmp = [];
        for i = 1: size(spks,2)
            tmp = [tmp,spks{cluster,i}]; 
        end

        [N,~] = histcounts(tmp,edges);
        res_sm(:,1) = movmean(N, mov_mean);

        for i = num_beats: size(spks,2)
            tmp = [];
            for j = -num_beats+1: 0
                tmp = [tmp,spks{cluster,i+j}];
            end
            [N,~] = histcounts(tmp,edges);
            res_sm(:,i-num_beats+2) = movmean(N, mov_mean);
            [~,idx] = max(res_sm(int(1):int(2),i-num_beats+1));
            peak(i-num_beats+1)= edges(int(1)+idx);
        
        end
        app.sldr_spikelag.Limits = [2 , size(res_sm,2)];
        app.sldr_spikelag.Value = 2;
        app.spike_res.spike_lag.res = res_sm;
        app.spike_res.spike_lag.peak = peak;

        app.sbtn_spklag.Value = true;
        app.sbtn_skew_kurt.Value = false;
        hb_ts = app.hb_res.t_events(app.hb_res.use_beats(:,int_idx));
        p=pcolor(app.ax_lag,hb_ts(num_beats:end), edges(2:end),res_sm(:,2:end));
        p.EdgeColor = 'interp';
        line (app.ax_lag,[hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
        line (app.ax_lag,hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
        app.spike_res.spike_lag.line =line (app.ax_lag,[-1, -1], [edges(2), edges(end)], 'Color', 'r','LineStyle',':');
        xlim (app.ax_lag,[hb_ts(num_beats),hb_ts(end)])

        plot (app.ax_hist, edges(2:end), app.spike_res.spike_lag.res(:,1))
        xlim (app.ax_hist, [edges(2),edges(end)])
        [~,i] = max(res_sm(int(1):int(2),1));
        tmp = edges(int(1):int(2));
        lag = mean([tmp(i),tmp(i+1)]);
        title(app.ax_hist,['global lag: ' num2str(lag) ' [s]'  ])
    
    case 'pointer'
        
        plot (app.ax_hist, edges(2:end), app.spike_res.spike_lag.res(:,round(app.sldr_spikelag.Value)))
        xlim (app.ax_hist, [edges(2),edges(end)])
        app.spike_res.spike_lag.line.XData = [app.hb_res.t_events(app.edt_num_beat.Value+round(app.sldr_spikelag.Value)-2),...
                                              app.hb_res.t_events(app.edt_num_beat.Value+round(app.sldr_spikelag.Value)-2)];
    case 'switch'


end














% tmp =[];
% for i = 1: length(spks)
%    tmp = [tmp,spks{i}]; 
% end
% edges = -0.5 :0.005:3;
% h =histogram(tmp,edges);
% 
% 
% 
% num_beats = 10;
% res = nan(length(edges)-1,length(spks)-num_beats+1);
% res_sm= nan(length(edges)-1,length(spks)-num_beats+1);
% peak = nan(1,length(spks)-num_beats+1);
% int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
% for i = num_beats: length(spks)
%     tmp = [];
%     for j = -num_beats+1: 0
%         tmp = [tmp,spks{i+j}];
%     end
%     [N,~] = histcounts(tmp,edges);
%     res(:,i-num_beats+1) = N;%h.BinCounts;
%     res_sm(:,i-num_beats+1) = movmean(res(:,i-num_beats+1), 20);
%     [~,idx] = max(res_sm(int(1):int(2),i-num_beats+1));
%     peak(i-num_beats+1)= edges(int(1)+idx);
% 
% end
% 
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res_sm);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)function update_spikelag_display(app, status)
% spks = app.spike_res.spike_lag.spks;
% 
% edges = -0.5 :0.005:3;
% 
% if strcmp(status, 'new')
%     num_beats = app.edt_num_beat.Value;
%         mov_mean = app.edt_lag_movmean.Value;
%         cluster = find(strcmp(app.popup_cluster_lag.Value,app.popup_cluster_lag.Items));
%         res_sm= nan(length(edges)-1,length(spks)-num_beats+2);
%         peak = nan(1,length(spks)-num_beats+1);
%         int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
%         tmp = [];
%         for i = 1: length(spks)
%             tmp = [tmp,spks{cluster,i}]; 
%         end
% 
%         [N,~] = histcounts(tmp,edges);
%         res_sm(:,1) = movmean(N, mov_mean);
% 
%         for i = num_beats: length(spks)
%             tmp = [];
%             for j = -num_beats+1: 0
%                 tmp = [tmp,spks{i+j}];
%             end
%             [N,~] = histcounts(tmp,edges);
%             res_sm(:,i-num_beats+2) = movmean(N, mov_mean);
%             [~,idx] = max(res_sm(int(1):int(2),i-num_beats+1));
%             peak(i-num_beats+1)= edges(int(1)+idx);
%         
%         end
%         app.sldr_spikelag.Limits = [2 , size(res_sm,2)];
%         app.sldr_spikelag.Value = 2;
% end
% hb_ts = app.hb_res.t_events;
% hist_idx = app.sldr_spikelag.Value;
% 
% if app.sbtn_spklag.Value 
%     p=pcolor(app.ax_lag,hb_ts(num_beats:end), edges(2:end),res_sm(:,2:end));
%     p.EdgeColor = 'interp';
%     line (app.ax_lag,[hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
%     line (app.ax_lag,hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
%     
%     plot (app.ax_hist, edge(2:end), res_sm(:,hist_idx))
% else
% 
% end
% 
% 
% 
% 
% bla = true 
% 
% 
% 
% 
% 
% tmp =[];
% for i = 1: length(spks)
%    tmp = [tmp,spks{i}]; 
% end
% edges = -0.5 :0.005:3;
% h =histogram(tmp,edges);
% 
% 
% 
% num_beats = 10;
% res = nan(length(edges)-1,length(spks)-num_beats+1);
% res_sm= nan(length(edges)-1,length(spks)-num_beats+1);
% peak = nan(1,length(spks)-num_beats+1);
% int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
% for i = num_beats: length(spks)
%     tmp = [];
%     for j = -num_beats+1: 0
%         tmp = [tmp,spks{i+j}];
%     end
%     [N,~] = histcounts(tmp,edges);
%     res(:,i-num_beats+1) = N;%h.BinCounts;
%     res_sm(:,i-num_beats+1) = movmean(res(:,i-num_beats+1), 20);
%     [~,idx] = max(res_sm(int(1):int(2),i-num_beats+1));
%     peak(i-num_beats+1)= edges(int(1)+idx);
% 
% end
% 
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res_sm);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
% 
% 
% num_beats = 10;
% res = nan(length(edges)-1,length(spks)-num_beats+1);
% res_sm= nan(length(edges)-1,length(spks)-num_beats+1);
% peak = nan(1,length(spks)-num_beats+1);
% int= [find(edges >= 0.7,1), find(edges >= 1.6,1) ];
% for i = num_beats: length(spks)
%     tmp = [];
%     for j = -num_beats+1: 0
%         tmp = [tmp,spks{i+j}];
%     end
%     [N,~] = histcounts(tmp,edges);
%     res(:,i-num_beats+1) = N;%h.BinCounts;
%     res_sm(:,i-num_beats+1) = movmean(res(:,i-num_beats+1), 20);
%     [~,idx] = max(res_sm(int(1):int(2),i-num_beats+1));
%     peak(i-num_beats+1)= edges(int(1)+idx);
% 
% end
% 
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)
% 
% figure
% p=pcolor(hb_ts(num_beats:end), edges(2:end),res_sm);
% p.EdgeColor = 'interp';
% line ([hb_ts(num_beats) hb_ts(end)], [0,0],'Color', 'r', 'LineWidth',2)
% line (hb_ts(num_beats:end), peak,'Color', 'g', 'LineWidth',2)