function update_transduction_axis(app)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

int_idx = find(strcmp(app.popup_int_transduction.Value,app.popup_int_transduction.Items));
int_name = app.popup_int_entropy.Value;
%sel_cluster = find(strcmp(app.popup_cluster.Value,app.popup_cluster.Items));

cla(app.ax_spike_shape)
cla(app.ax_trans_1_1)
cla(app.ax_trans_2_1)
cla(app.ax_trans_1_2)
cla(app.ax_trans_2_2)
cla(app.ax_trans_1_3)
cla(app.ax_trans_2_3)
leg_str = {[]};
yinf=-1;
ysup=1;
for i=1 :length(app.spike_res.transduction(int_idx).cluster)
    hold(app.ax_spike_shape,'on')
    plot (app.ax_spike_shape, app.spike_res.transduction(int_idx).cluster(i).mean_shape )
    hold(app.ax_spike_shape,'on')
    leg_str{i} = ['cl: ' num2str(i) ' n=' num2str(app.spike_res.transduction(int_idx).cluster(i).n_spikes)];

    x = [-0.5 0.5 0.5 -0.5];
    y = [-1 -1 1 1];
    hold(app.ax_trans_1_1, 'on')
    plot(app.ax_trans_1_1, app.spike_res.transduction(int_idx).cluster(i).fr_bp(:,1), app.spike_res.transduction(int_idx).cluster(i).fr_bp(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson FR - BP'), xlabel('cardiac cycle'),title('Pearson FR - BP')
    patch(app.ax_trans_1_1, x, y, [0.7 0.7 0.7])
    app.ax_trans_1_1.Toolbar = [];
    ylim(app.ax_trans_1_1, [yinf ysup])
    hold(app.ax_trans_1_1, 'off')
    
    hold(app.ax_trans_2_1, 'on')
    plot(app.ax_trans_2_1, app.spike_res.transduction(int_idx).cluster(i).fr_rr(:,1), app.spike_res.transduction(int_idx).cluster(i).fr_rr(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson FR - RR'), xlabel('cardiac cycle'),title('Pearson FR - RR')
    patch(app.ax_trans_2_1, x, y, [0.7 0.7 0.7])
    app.ax_trans_2_1.Toolbar = [];
    ylim(app.ax_trans_2_1,[yinf ysup])
    hold(app.ax_trans_2_1, 'off')
    
    hold(app.ax_trans_1_2, 'on')
    plot(app.ax_trans_1_2, app.spike_res.transduction(int_idx).cluster(i).latency_bp(:,1), app.spike_res.transduction(int_idx).cluster(i).latency_bp(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson First Spike Lat - BP'), xlabel('cardiac cycle'),title('Pearson First Spike Lat - BP')
    patch(app.ax_trans_1_2, x, y, [0.7 0.7 0.7])
    app.ax_trans_1_2.Toolbar = [];
    ylim(app.ax_trans_1_2, [yinf ysup])
    hold(app.ax_trans_1_2, 'off')
    
    hold(app.ax_trans_2_2, 'on')
    plot(app.ax_trans_2_2, app.spike_res.transduction(int_idx).cluster(i).latency_rr(:,1), app.spike_res.transduction(int_idx).cluster(i).latency_rr(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson First Spike Lat - RR'), xlabel('cardiac cycle')
    patch(app.ax_trans_2_2, x, y, [0.7 0.7 0.7])
    app.ax_trans_2_2.Toolbar = [];
    ylim(app.ax_trans_2_2, [yinf ysup])
    hold(app.ax_trans_2_2, 'off')
    
    hold(app.ax_trans_1_3, 'on')
    plot(app.ax_trans_1_3, app.spike_res.transduction(int_idx).cluster(i).phase_bp(:,1), app.spike_res.transduction(int_idx).cluster(i).phase_bp(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson First Spike PHASE - BP'), xlabel('cardiac cycle'),title('Pearson First Spike PHASE - BP')
    app.ax_trans_1_3.Toolbar = [];
    ylim(app.ax_trans_1_3, [yinf ysup])
    hold(app.ax_trans_1_3, 'off')

    hold(app.ax_trans_2_3, 'on')
    plot(app.ax_trans_2_3, app.spike_res.transduction(int_idx).cluster(i).phase_rr(:,1), app.spike_res.transduction(int_idx).cluster(i).phase_rr(:,2),'-o','linewidth',1.5)
%     ylabel('Pearson First Spike PHASE - RR'), xlabel('cardiac cycle'),title('Pearson First Spike PHASE - RR')
    app.ax_trans_2_3.Toolbar = [];
    ylim(app.ax_trans_2_3, [yinf ysup])
    hold(app.ax_trans_2_3, 'off')
end
lgnd = legend(app.ax_spike_shape,leg_str,'Location','southeast');
lgnd.ItemTokenSize(1) = 10;
end

