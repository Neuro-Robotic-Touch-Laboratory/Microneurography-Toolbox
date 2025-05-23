function set_burst_anal(app)
%SET_BURST_ANAL checks if heatbeats / respiration data is available and
%sets the burst analysis type
if  app.sbtn_thresh.Value
    if app.sbtn_shape.Value
        type = false;
    else
        type = true;
    end
else
    type = false;
end

if  app.sbtn_resp.Value
    if app.sbtn_hb.Value
        type(2) = false;
    else
        type(2) = true;
    end
else
    type(2) = false;
end

if isempty(app.hb_res)
    hb = false;
else
    hb = true;
end
if isempty(app.resp_res)
    resp = false;
else
    resp = true;
end

if ~hb & ~resp 
    app.pnl_burst_shape.Visible = "off";
    app.pnl_burst_shape.Enable = "off";
    app.pnl_burst_thresh.Visible = "on";
    app.pnl_burst_thresh.Enable = "on";
    app.sbtn_thresh.Value = true;
    app.sbtn_shape.Value = false;
    app.sbtn_shape.Enable = "off";
else
    app.sbtn_shape.Enable = "on";
    if type(1)
        app.pnl_burst_shape.Visible = "off";
        app.pnl_burst_shape.Enable = "off";
        app.pnl_burst_thresh.Visible = "on";
        app.pnl_burst_thresh.Enable = "on";
        app.sbtn_thresh.Value = true;
        app.sbtn_shape.Value = false;
        

    else
        app.pnl_burst_shape.Visible = "on";
        app.pnl_burst_shape.Enable = "on";
        app.pnl_burst_thresh.Visible = "off";
        app.pnl_burst_thresh.Enable = "off";
        app.sbtn_thresh.Value = false;
        app.sbtn_shape.Value = true;
        
        if resp & hb
            app.sbtn_resp.Enable = "on";
            app.sbtn_hb.Enable = "on";
            if type(2)
                app.sbtn_resp.Value = true;
                app.sbtn_hb.Value = false;
            else
                app.sbtn_resp.Value = false;
                app.sbtn_hb.Value = true;
            end
        else
            if hb 
                app.sbtn_hb.Enable = "on";
                app.sbtn_resp.Enable = "off";
                app.sbtn_resp.Value = false;
                app.sbtn_hb.Value = true;
            else
                app.sbtn_hb.Enable = "off";
                app.sbtn_resp.Enable = "on";
                app.sbtn_resp.Value = true;
                app.sbtn_hb.Value = false;
            end
        end

    
    end
end



end

