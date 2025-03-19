function [RRR7, RRR8 ] = transduction_calc_v2(data, ts_ecg_in, ts_bpfoot_in, bpValues, t_spikes, idx_bpfoot_in)


ts_ecg = ts_ecg_in;
ts_bpfoot = ts_bpfoot_in;
idx_bpfoot = idx_bpfoot_in;

tmp = find (ts_bpfoot < ts_ecg(1));

ts_bpfoot(tmp) = [];
idx_bpfoot(tmp) = [];

[ts_ecg_cut, ts_bp_cut, idx_bp_cut] =  realign_bp_ecg(ts_ecg, ts_bpfoot, idx_bpfoot);







end
 