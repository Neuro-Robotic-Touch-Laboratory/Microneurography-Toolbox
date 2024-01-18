function auto_resp(app)
%auto_resp checks if respiration is already detected and runs  if
%not 
%   Detailed explanation goes here

if isempty(app.resp_res)
    
    [data,ts,~, ~] = current_signal(app, app.settings.channel_idx.resp);
    data = lowpass(data,3,1/ts(1));
    data(:,2) = (ts(1):ts(1):ts(2))';
    app.resp_res = resp_det2(app,data, .25);
    app.btn_resp.BackgroundColor = 'g';
    derived_resp_signals(app)
    app.settings.done.resp = true;

end
end