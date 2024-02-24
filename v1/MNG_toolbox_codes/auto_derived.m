function auto_derived(app)
%AUTO_DERIVED Summary of this function goes here
%   Detailed explanation goes here
if ~isempty(app.data)
    if isempty(app.burst_ints)
        [~,  idx] = max(vertcat(app.data.tic_multipl));
        [~,ts,~, ~] = current_signal(app, idx);
        bords = [ts(1), ts(2)-ts(1)];
        app.burst_ints = struct('name','autofull','type', 1, 'parent', 'none', 'borders', bords);
        app.btn_sel_int.BackgroundColor = [0,1,0];
    end
    
    auto_bp(app)
    auto_hb(app)
    auto_resp(app)
    auto_etco2(app)
%     auto_min_vent(app)
else
    warndlg('Please load datafile first')
end

end

