function auto_etco2(app)
%AUTO_ETCO2 Summary of this function goes here
%   Detailed explanation goes here
if ~isnan(app.settings.channel_idx.co2)
    if isempty(find(contains(vertcat(app.data.name),'etCO2')))
        [data,ts,~, ~] = current_signal(app, app.settings.channel_idx.co2);
        t= ts(1):ts(1):ts(2);
        [b,a] =  butter (2, 10/(1/ts(1)));
        data_flt = filtfilt(b,a,data);
        [pks,lcs] = findpeaks(data_flt,'MinPeakDistance',1000,'MinPeakProminence',std(data_flt));
        pks(2:end+1) = pks;
        pks(end+1)= pks(end);
        lcs = [1;lcs;length(data_flt)];
        data_co2 = interp1(t(lcs),pks,t);
        app.data(end+1).name = 'etCO2';
        app.data(end).unit = 'mmHg';
        app.data(end).tic_multipl = app.data(app.settings.channel_idx.co2).tic_multipl;
        app.data(end).data = data_co2';
        app.data(end).derived = true;
        app.data(end).ts = ts;
    end 
end
end

