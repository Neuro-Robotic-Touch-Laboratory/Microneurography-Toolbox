classdef pointer_resp < handle
    %POINTER manages  the data to display and the rwaves to display
    %   Detailed explanation goes here
    
    properties
        start_idx
        stop_idx
        dspl_strt_idx
        dspl_peak_idx
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
        function obj = pointer_resp (app)
            obj.dur = int64(app.edt_window.Value*1/(mean(diff(app.data(:,2))))); 
            obj.start_idx = int64(1);
            obj.stop_idx = int64(obj.start_idx +obj.dur);
            obj.dspl_strt_idx = app.settings.start_idx(app.settings.start_idx >= obj.start_idx & app.settings.start_idx <= obj.stop_idx);
            obj.dspl_peak_idx = app.settings.peak_idx(app.settings.peak_idx >= obj.start_idx & app.settings.peak_idx <= obj.stop_idx);
            obj.ax_1 = app.ax_ov_resp;
            obj.ax_2 = app.ax_edit;
            obj.length = size(app.data,1);
            notify (obj, 'new_frame')
            
        end
        
        function  obj = next_frame(obj,app)
            if obj.stop_idx -1000+obj.dur > obj.length
                obj.stop_idx = int64(obj.length);
                obj.start_idx = int64(obj.stop_idx-obj.dur);
            else 
                obj.start_idx = int64(obj.stop_idx-1000);
                obj.stop_idx = int64(obj.start_idx +obj.dur);
            end
            obj.dspl_strt_idx = app.settings.start_idx(app.settings.start_idx >= obj.start_idx & app.settings.start_idx <= obj.stop_idx);
            obj.dspl_peak_idx = app.settings.peak_idx(app.settings.peak_idx >= obj.start_idx & app.settings.peak_idx <= obj.stop_idx);
            
            notify (obj,'new_frame')
        end
        
         function  obj = prev_frame(obj,app)
            if obj.start_idx +1000-obj.dur < 1
                obj.start_idx = int64(1);
                obj.stop_idx = int64(obj.start_idx +obj.dur);
            else 
                obj.start_idx = int64(obj.start_idx+1000 -obj.dur);
                obj.stop_idx = int64(obj.start_idx +obj.dur);
            end
            obj.dspl_strt_idx = app.settings.start_idx(app.settings.start_idx >= obj.start_idx & app.settings.start_idx <= obj.stop_idx);
            obj.dspl_peak_idx = app.settings.peak_idx(app.settings.peak_idx >= obj.start_idx & app.settings.peak_idx <= obj.stop_idx);
            notify (obj,'new_frame')
        end
               
        function obj = get_frame(obj, idx,app)
            obj.start_idx = int64(idx);
            obj.stop_idx = obj.start_idx +obj.dur;
            if obj.start_idx < 1
                obj.start_idx = int64(1);
                obj.stop_idx = int64(obj.start_idx +obj.dur);
            end
            if obj.stop_idx > obj.length
                obj.stop_idx = int64(obj.length);
                obj.start_idx = int64(obj.stop_idx-obj.dur);
            end
            obj.dspl_strt_idx = app.settings.start_idx(app.settings.start_idx >= obj.start_idx & app.settings.start_idx <= obj.stop_idx);
            obj.dspl_peak_idx = app.settings.peak_idx(app.settings.peak_idx >= obj.start_idx & app.settings.peak_idx <= obj.stop_idx);
            notify(obj, 'new_frame')
        end
            
    end
    
end

