function auto_bp(app) 
%auto_Bb checks if bb is already detected and runs  if
%not 
%   Detailed explanation goes here
if isempty(app.bp_res)
    if ~isnan(app.settings.channel_idx.bldp)
        [data,ts,~, ~] = current_signal(app, app.settings.channel_idx.bldp);
        data = lowpass(data,3,1/ts(1));
        try
            [footIndex, systolicIndex, notchIndex, dicroticIndex ] = bp_dect( data, 200, 1, 'mmHg', 1);
        catch
            warndlg('no proper bloodpressure waveforms detected - bloodpressure calculation skipped')
            app.bp_res.fail = true;
            return
        end
        results.foot_idx = unique(footIndex);
        results.systolic_idx = unique(systolicIndex);
        results.notch_idx = unique(notchIndex);
        results.dicrotic_idx = unique(dicroticIndex);
        results.ts = ts(1);
        app.bp_res = results;
        derived_bp_signals(app)
    end
end
end