function auto_hb(app)
%auto_hb checks if hb is already detected and runs hb_det_2 on the ecg if
%not 
%   Detailed explanation goes here

if isempty(app.hb_res)
    if ~isnan(app.settings.channel_idx.ecg)
        if app.settings.interval(1,1) < app.data(app.settings.channel_idx.ecg).ts(1)
            app.settings.interval(1,1) = app.data(app.settings.channel_idx.ecg).ts(1);
        end
        if app.settings.interval(1,2) > app.data(app.settings.channel_idx.ecg).ts(1)*length(app.data(app.settings.channel_idx.ecg).data)
            app.settings.interval(1,2) = app.data(app.settings.channel_idx.ecg).ts(1)*length(app.data(app.settings.channel_idx.ecg).data);
        end
    
        ecg_r = app.data(app.settings.channel_idx.ecg).data;
        ecg_use = ecg_r(app.settings.interval(1,1)/app.data(app.settings.channel_idx.ecg).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.ecg).ts(1));
        ecg(:,1) = lowpass(ecg_use(:,1),40,1/app.data(app.settings.channel_idx.ecg).ts(1));
        ecg(:,2) = (0:size(ecg_use,1)-1)*app.data(app.settings.channel_idx.ecg).ts(1);
        if isempty(app.bp_res)
            bp = [];
            use_bp = false;
        else
            bp = app.bp_res;
            use_bp = true;
            if isfield(bp ,'fail')
                if bp.fail
                    use_bp = false;
                end
            end
            
        end
        
        app.hb_res = hb_det_2(ecg,3,false,bp,use_bp);
        app.hb_res.dt_instantaneous = diff(app.hb_res.t_events);
        app.hb_res.rate_instantaneous = (1/app.hb_res.dt_instantaneous * 60)';
        app.hb_res.rate_mean = mean(app.hb_res.rate_instantaneous);
        app.hb_res.dtDifferences = diff(app.hb_res.dt_instantaneous);
        app.hb_res.ts = [ecg(1,2),mean(diff(ecg(:,2))), ecg(end,2)];
        app.btn_hb.BackgroundColor = 'g';
    
        derived_cardiac_signals(app)
        app.settings.done.hb = true;
       
    end
    sel_beats(app)
end

end
