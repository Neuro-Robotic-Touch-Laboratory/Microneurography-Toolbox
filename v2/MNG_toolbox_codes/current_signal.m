function [data,ts,name, unit] = current_signal(app, channel_index)
%CURRENT_SIGNAL returns data, dt and max t , name and units for the current
%interval for the selected signal. adapts for the different timescale of
%derived signals
%   Detailed explanation goes here

if app.data(channel_index).derived
    data = app.data(channel_index).data;
else
    
%     data = app.data(channel_index).data(int32(app.settings.interval(1,1)/app.data(channel_index).ts(1) +1): ...
%                                         int32(app.settings.interval(1,2)/app.data(channel_index).ts(1)));

    strt = int32(app.settings.interval(1,1)/app.data(channel_index).ts(1) +1);
    stp = int32(app.settings.interval(1,2)/app.data(channel_index).ts(1));
    if strt <1
        strt = 1;
    end
    if stp > length(app.data(channel_index).data)
        stp = length(app.data(channel_index).data);
    end
    data = app.data(channel_index).data(strt: stp);
end

ts = [1, length(data)]*app.data(channel_index).ts(1);
name = app.data(channel_index).name;
unit = app.data(channel_index).unit;
end

