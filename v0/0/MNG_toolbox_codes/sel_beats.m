function sel_beats(app)
%SEL_BEATS Summary of this function goes here
%   Detailed explanation goes here
if  ~isempty(app.burst_ints) && ~isempty(app.hb_res)
    app.hb_res.use_beats = [];
    ts_beats = app.hb_res.ts(1) : app.hb_res.ts(2) : app.hb_res.ts(3);
    app.hb_res.use_beats = true(size(app.hb_res.idx));
    
    for i = 1: length(app.burst_ints)
    
        app.hb_res.use_beats(1:length(app.hb_res.idx),i+1) = false(size(app.hb_res.idx));
        for j = 1: size(app.burst_ints(i).borders)
            ts_idx = [find(ts_beats >= app.burst_ints(i).borders(j,1),1),...
                      find(ts_beats <= app.burst_ints(i).borders(j,2),1,'last')];
            tmp_idx(1) = find(app.hb_res.idx >= ts_idx(1),1);
            tmp_idx(2) = find(app.hb_res.idx <= ts_idx(2),1,'last');
            app.hb_res.use_beats(tmp_idx(1):tmp_idx(2),i+1) = true;
        end
    
    end
end
end

