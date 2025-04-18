function auto_min_vent(app)

if ~isnan(app.settings.channel_idx.res_flow)
    if isempty(find(contains(vertcat(app.data.name),'min. vent.')))
        [data,ts,~, ~] = current_signal(app, app.settings.channel_idx.res_flow);
        t= ts(1):ts(1):ts(2);
        sig_p = data;
        sig_p(sig_p<0) = 0;
        sig_n = data;
        sig_n(sig_p>0) = 0;
        
        win = 60*(round(1/ts(1)));
        ti = t(round(win /2)+ (1 : round((1/ts(1))/5) : length(data)-win));
        int_p = nan(size(ts));
        int_n = nan(size(ts));
        idx = 1 : round((1/ts(1))/5) : length(data)-win;
        for i = 1:length(idx)
            int_p(i) = trapz(sig_p(idx(i):idx(i)+win -1));
            int_n(i) = trapz(abs(sig_n(idx(i):idx(i)+win -1)));
        end
        rat = int_p./int_n;
        int_p(rat>= 1.2) = 0;
        int_p(rat<= 0.8) = 0;
        int_p = [int_p(1) ,int_p, int_p(end)];
        ti = [t(1), ti,t(end)];
        int_p(int_p ==0 ) = nan;
        min_vent = interp1(ti,int_p,t);

        app.data(end+1).name = 'min. vent.';
        app.data(end).unit = 'l/min';
        app.data(end).tic_multipl = app.data(app.settings.channel_idx.res_flow).tic_multipl;
        app.data(end).data = min_vent;
        app.data(end).derived = true;
        app.data(end).ts = ts;
       
    end 
end
end

