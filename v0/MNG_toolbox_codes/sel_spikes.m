function spike_res_out = sel_spikes(spike_res_in, ts, burst_ints)
%SEL_SPIKES Summary of this function goes here
%   Detailed explanation goes here


spike_res_in.use_spikes = [];
ts_msna  = ts(1) : ts(1) : ts(2);
spike_res_in.use_spikes(1:length(spike_res_in.spike_idx),1) = ones(size(spike_res_in.spike_idx));
for i =1 : length(burst_ints)
    spike_res_in.use_spikes(1:length(spike_res_in.spike_idx),i+1) = zeros(size(spike_res_in.spike_idx));
    for j = 1: size(burst_ints(i).borders)
        ts_idx = [find(ts_msna >= burst_ints(i).borders(j,1),1),...
                  find(ts_msna <= burst_ints(i).borders(j,2),1,'last')];
        tmp_idx(1) = find(spike_res_in.spike_idx >= ts_idx(1),1);
        tmp_idx(2) = find(spike_res_in.spike_idx <= ts_idx(2),1,'last');
        spike_res_in.use_spikes(tmp_idx(1):tmp_idx(2),i+1) = 1;
    end
end

spike_res_out = spike_res_in;
end

