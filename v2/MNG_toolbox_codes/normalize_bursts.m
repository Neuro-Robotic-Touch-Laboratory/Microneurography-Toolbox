function  normalize_bursts(app)
%NORMALIZE_BURSTS Summary of this function goes here
%   normalizes the burst amplitude ant integral by the amplitude 
%   of the highest burst in the "baseline_interval" 
ts = app.burst_res.ts(1):app.burst_res.ts(2):app.burst_res.ts(3);
if app.chkbx_normalize.Value

    bsln_int = find(strcmp(app.popup_baseline_int.Value, app.popup_baseline_int.Items))+1;
    
    bl_bursts_idx = find(app.burst_res.use_burst(:,1) & app.burst_res.use_burst(:,bsln_int));
    
    bl_burst_amp =nan(1,length(bl_bursts_idx)); 
    
    
    for i = 1:length(bl_bursts_idx)
        switch app.burst_res.type
            case 1
                tmp_d = [ts(app.burst_res.burst_loc(bl_bursts_idx(i),1):app.burst_res.burst_loc(bl_bursts_idx(i),2))',app.burst_res.x(app.burst_res.burst_loc(bl_bursts_idx(i),1):app.burst_res.burst_loc(bl_bursts_idx(i),2))'*1000];      
                tmp_d(:,2) = tmp_d(:,2) - min (tmp_d(:,2));
                
                bl_burst_amp(i) = max(tmp_d(:,2));
            case 2
                bl_burst_amp(i) = ampl_tresh(app.burst_res.x, app.burst_res.burst_loc(bl_bursts_idx(i),:), 1);
        end
    end
    
    factor = 100/max(bl_burst_amp);

else
    factor = 1;

end

app.burst_res.burst_int = nan(size(app.burst_res.burst_int));
app.burst_res.burst_amp = nan(size(app.burst_res.burst_amp));

for i =1:size(app.burst_res.burst_loc,1)
    switch app.burst_res.type
        case 1
            tmp_d = [ts(app.burst_res.burst_loc(i,1):app.burst_res.burst_loc(i,2))',app.burst_res.x(app.burst_res.burst_loc(i,1):app.burst_res.burst_loc(i,2))'*1000];
            tmp_d(:,2) = tmp_d(:,2) *factor;
            tmp_d(:,2) = tmp_d(:,2) - min (tmp_d(:,2));
            
            app.burst_res.burst_amp(i) = max(tmp_d(:,2));
            app.burst_res.burst_int(i) = trapz(tmp_d(:,1),tmp_d(:,2));
        case 2
            app.burst_res.burst_amp(i) = ampl_tresh(app.burst_res.x, app.burst_res.burst_loc(i,:), factor);
            app.burst_res.burst_int(i) = integal_tresh(app.burst_res.xshift,ts, app.burst_res.burst_loc(i,:), factor);
    end
end



end

function amp = ampl_tresh(data, loc, factor)
    amp = max(data(loc(1):loc(2)-50))*factor;
end


function integral = integal_tresh(data, ts, loc, factor)
    integral = trapz(ts(loc(1):loc(2)),data(loc(1):loc(2))*factor);
end

