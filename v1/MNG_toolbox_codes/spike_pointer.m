classdef spike_pointer < handle
    %POINTER manages  the data to display and the rwaves to display
    %   Detailed explanation goes here
    
    properties
        start_idx_det
        stop_idx_det
        start_idx_ov
        stop_idx_ov
        dur_det
        dur_ov
        length
        overlap
        plot_idx
    end
    
    events
        new_frame
        new_frame_ov
    end
    
    methods
        function obj = spike_pointer (app)
            obj.dur_det = int64(app.edt_framesize_det.Value*1/(mean(diff(app.data(:,2))))); 
            obj.overlap = int64(round(obj.dur_det*0.15));
            obj.dur_ov = int64(app.edt_framesize_ov.Value*1/(mean(diff(app.data(:,2)))));
            
            obj.start_idx_det = int64(1);
            obj.stop_idx_det = obj.start_idx_det+obj.dur_det;
            obj.start_idx_ov = int64(1);
            obj.stop_idx_ov =  obj.start_idx_ov + obj.dur_ov;
            obj.length = size(app.data,1);
%             spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
%             displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
%             obj.plot_idx = {[]};
%             for i =1:app.spike_res.n_clus+1
%                 if i ==1
%                     obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
%                 else
%                     obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
%                 end
%             end
%             notify (obj, 'new_frame_ov')
            
        end

        function obj = set_ov_dur(obj,dur,app)
            obj.dur_ov = int64(dur*1/mean(diff(app.data(:,2))));
            obj.stop_idx_ov = obj.start_idx_ov+obj.dur_ov;
            if obj.start_idx_ov < 1
                obj.start_idx_ov = int64(1);
                obj.stop_idx_ov = int64(obj.start_idx_ov +obj.dur_ov);
            end
            if obj.stop_idx_ov > obj.length
                obj.stop_idx_ov = int64(obj.length);
                obj.start_idx_ov = int64(obj.stop_idx_ov -obj.dur_ov);
            end
            obj.start_idx_det = obj.start_idx_ov;
            obj.stop_idx_det = obj.start_idx_det+obj.dur_det;
            app.sldr_full.Enable = 'on';
            app.sldr_full.Limits = [double(1),double(size(app.data,1)-app.ptr.dur_ov)];
            app.sldr_full.Value = obj.start_idx_ov;
            app.sldr_overview.Limits = [double(obj.start_idx_ov), double(obj.stop_idx_ov-obj.dur_det)];
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            notify(obj, 'new_frame_ov')
        end
        
        function obj = set_det_dur(obj,dur,app)
            obj.dur_det = int64(dur*1/mean(diff(app.data(:,2))));
            obj.stop_idx_det = obj.start_idx_det +obj.dur_det;
            if obj.stop_idx_det > obj.stop_idx_ov
                obj.stop_idx_det = obj.stop_idx_ov;
                obj.start_idx_det = int64(obj.stop_idx_det -obj.dur_det);
            end
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            app.sldr_overview.Limits = [double(obj.start_idx_ov), double(obj.stop_idx_ov-obj.dur_det)];
            
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            notify(obj, 'new_frame')
        end

        function  obj = next_frame(obj,app)
            if obj.stop_idx_det -obj.overlap +obj.dur_det > obj.stop_idx_ov
                obj.stop_idx_det = obj.stop_idx_ov;
                obj.start_idx_det = int64(obj.stop_idx_det-obj.dur_det);
            else 
                obj.start_idx_det = int64(obj.stop_idx_det-obj.overlap);
                obj.stop_idx_det = int64(obj.start_idx_det +obj.dur_det);
            end
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            notify (obj,'new_frame')
        end
        
         function  obj = prev_frame(obj,app)
            if obj.start_idx_det +obj.overlap-obj.dur_det < obj.start_idx_ov
                obj.start_idx_det = obj.start_idx_ov;
                obj.stop_idx_det = int64(obj.start_idx_det +obj.dur_det);
            else 
                obj.start_idx_det = int64(obj.start_idx_det +obj.overlap -obj.dur_det);
                obj.stop_idx_det = int64(obj.start_idx_det +obj.dur_det);
            end
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            notify (obj,'new_frame')
        end
               
        function obj = get_frame(obj, idx,app)
            obj.start_idx_det = int64(idx);
            obj.stop_idx_det = obj.start_idx_det +obj.dur_det;
            if obj.start_idx_det < obj.start_idx_ov
                obj.start_idx_det = obj.start_idx_ov;
                obj.stop_idx_det = int64(obj.start_idx_det +obj.dur_det);
            end
            if obj.stop_idx_det > obj.stop_idx_ov
                obj.stop_idx_det = obj.stop_idx_ov;
                obj.start_idx_det = int64(obj.stop_idx_det-obj.dur_det);
            end
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            notify(obj, 'new_frame')
        end

        function obj = get_frame_ov(obj, idx,app)
            obj.start_idx_ov = int64(idx);
            obj.stop_idx_ov = idx +obj.dur_ov;
            if obj.start_idx_ov < 1
                obj.start_idx_ov = int64(1);
                obj.stop_idx_ov = int64(obj.start_idx_ov +obj.dur_ov);
            end
            if obj.stop_idx_ov > obj.length
                obj.stop_idx_ov = int64(obj.length);
                obj.start_idx_ov = int64(obj.stop_idx_ov -obj.dur_ov);
            end
            obj.start_idx_det = obj.start_idx_ov;
            obj.stop_idx_det = obj.start_idx_det+obj.dur_det;
            app.sldr_overview.Limits = [double(obj.start_idx_ov), double(obj.stop_idx_ov-obj.dur_det)];
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_det) & (app.spike_res.spike_idx<= obj.stop_idx_det));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.plot_idx = {[]};
            for i =1:app.spike_res.n_clus+1
                if i ==1
                    obj.plot_idx{i} = displ_spks(~displ_spks(:,2),1);
                else
                    obj.plot_idx{i} = displ_spks(displ_spks(:,2) & displ_spks(:,3) == i-1,1);
                end
            end
            
            notify(obj, 'new_frame_ov')
        end
            
    end
    
end

