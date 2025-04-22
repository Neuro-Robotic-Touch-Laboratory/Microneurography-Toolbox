function beat_spk_lag = get_spike_lag(spike_ts, spike_clust,spike_use, hb_ts, beats_use, int_idx)

spike_ts = spike_ts/1000;
num_clus = max(spike_clust);
spks = {};
spike_use = spike_use(:,int_idx);
hb_ts = hb_ts(beats_use(:,int_idx));
for i = 1 : length(hb_ts) 
    spks{1,i} = spike_ts(spike_ts>=(hb_ts(i)-0.5) & spike_ts<=(hb_ts(i)+3) & spike_use');
    spks{1,i} = spks{1,i}-hb_ts(i);
    for j= 1:num_clus
        spks{1+j,i} = spike_ts(spike_ts>=(hb_ts(i)-0.5) & spike_ts<=(hb_ts(i)+3) & spike_clust == j  & spike_use');
        spks{1+j,i} = spks{1+j,i}-hb_ts(i);
    end
end

beat_spk_lag = spks;