function [ts_ecg_out, ts_bp_out, idx_bp_out] =  realign_bp_ecg(ts_ecg_in, ts_bp_in, idx_bp_in)

ts_ecg = unique(ts_ecg_in);
ts_bp = unique(ts_bp_in);
idx_bp = unique(idx_bp_in);

tmp = find (diff(ts_bp)< 0.5 );
ts_bp(tmp) = [];
idx_bp(tmp) = [];

tmp = find(diff(ts_ecg) < 0.5);
ts_ecg(tmp) = [];
