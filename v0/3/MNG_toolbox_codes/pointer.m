classdef pointer < handle
    %POINTER manages  the data to display and the rwaves to display
    %   Detailed explanation goes here
    
    properties
        start_idx
        stop_idx
        dspl_idx
        pre_idx
        post_idx
        dur
        ax_1
        ax_2
        length
    end
    
    events
        new_frame
    end
    
    methods
        function obj = pointer (app)
            obj.dur = app.edt_window.Value*1/(mean(diff(app.ecg(:,2)))); 
            if obj.dur > size(app.ecg,1)
                obj.dur = size(app.ecg,1)-1;
                app.edt_window.Value = floor(obj.dur*(mean(diff(app.ecg(:,2)))));
            end
            obj.start_idx = 1;
            obj.stop_idx = obj.start_idx +obj.dur;
            obj.dspl_idx = app.hb_res.idx(app.hb_res.idx >= obj.start_idx & app.hb_res.idx <= obj.stop_idx);
            obj.ax_1 = app.ax_ov_ecg;
            obj.ax_2 = app.ax_edit;
            obj.length = size(app.ecg,1);
            notify (obj, 'new_frame')
            
        end
        
        function  obj = next_frame(obj,app)
            if obj.stop_idx -1000+obj.dur > obj.length
                obj.stop_idx = obj.length;
                obj.start_idx = obj.stop_idx-obj.dur;
            else 
                obj.start_idx = obj.stop_idx-1000;
                obj.stop_idx = obj.start_idx +obj.dur;
            end
            obj.start_idx = int64(obj.start_idx);
            obj.stop_idx = int64(obj.stop_idx);
            obj.dspl_idx = app.hb_res.idx(app.hb_res.idx >= obj.start_idx & app.hb_res.idx <= obj.stop_idx);
            notify (obj,'new_frame')
        end
        
         function  obj = prev_frame(obj,app)
            if obj.start_idx +1000-obj.dur < 1
                obj.start_idx = 1;
                obj.stop_idx = obj.start_idx +obj.dur;
            else 
                obj.start_idx = obj.start_idx+1000 -obj.dur;
                obj.stop_idx = obj.start_idx +obj.dur;
            end
            obj.start_idx = int64(obj.start_idx);
            obj.stop_idx = int64(obj.stop_idx);
            obj.dspl_idx = app.hb_res.idx(app.hb_res.idx >= obj.start_idx & app.hb_res.idx <= obj.stop_idx);
            notify (obj,'new_frame')
        end
               
        function obj = get_frame(obj, idx,app)
            obj.start_idx = idx;
            obj.stop_idx = obj.start_idx +obj.dur;
            if obj.start_idx  < 1
                obj.start_idx = 1;
                obj.stop_idx = obj.start_idx +obj.dur;
            end
            if  obj.stop_idx > obj.length
                obj.start_idx = obj.length-obj.dur;
                obj.stop_idx = obj.length;
            end
            obj.start_idx = int64(obj.start_idx);
            obj.stop_idx = int64(obj.stop_idx);
            obj.dspl_idx = app.hb_res.idx(app.hb_res.idx >= obj.start_idx & app.hb_res.idx <= obj.stop_idx);
            notify(obj, 'new_frame')
        end
            
    end
    
end

