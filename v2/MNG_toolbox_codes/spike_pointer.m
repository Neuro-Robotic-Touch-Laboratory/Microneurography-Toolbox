classdef spike_pointer < handle
    %POINTER manages  the data to display and the rwaves to display
    %   Detailed explanation goes here
    
    properties
        start_idx_det
        stop_idx_det
        start_idx_ov
        stop_idx_ov
        start_idx_full
        stop_idx_full
        dur_det
        dur_ov
        overlap
        plot_idx
        ov_idx %% set
        data_len
    end
    
    events
        new_frame
        new_frame_ov
        new_frame_full
    end
    
    methods
        function obj = spike_pointer (app,full_idx)
            obj.data_len = full_idx(2);
            obj.dur_det = int64(app.edt_framesize_det.Value*1/(mean(diff(app.data(:,2))))); 
            obj.overlap = int64(round(obj.dur_det*0.15));
            obj.dur_ov = int64(app.edt_framesize_ov.Value*1/(mean(diff(app.data(:,2)))));
            obj.start_idx_full = full_idx(1);
            obj.stop_idx_full = full_idx(2);
            if (full_idx(2)-full_idx(1)+1)< obj.dur_ov
                obj.dur_ov = full_idx(2)-full_idx(1)+1;
                app.edt_framesize_ov.Value = obj.dur_ov+mean(diff(app.data(:,2)));
            end
        
            obj.start_idx_det = int64(full_idx(1));
            obj.stop_idx_det = obj.start_idx_det+obj.dur_det;
            obj.start_idx_ov = int64(full_idx(1));
            obj.stop_idx_ov =  obj.start_idx_ov + obj.dur_ov;

            %% set ov_idx
            
            %obj.length = size(app.data,1);
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

        function obj = set_full_idx(obj, full_idx,app)
            obj.start_idx_full = full_idx(1);
            obj.stop_idx_full = full_idx(2);
            obj.start_idx_det = int64(full_idx(1));
            obj.stop_idx_det = obj.start_idx_det+obj.dur_det;
            if (obj.stop_idx_full-obj.start_idx_full+1)< obj.dur_ov
                obj.dur_ov = obj.stop_idx_full-obj.start_idx_full+1;
                app.edt_framesize_ov.Value = fix((obj.dur_ov*mean(diff(app.data(:,2))))*1000)/1000;
            end
            obj.start_idx_ov = int64(full_idx(1));
            obj.stop_idx_ov =  obj.start_idx_ov + obj.dur_ov;
            tmp = [full_idx(1), full_idx(2)-obj.dur_ov+1];
            if tmp(1) == tmp(2)
                app.sldr_full.Enable = 'off';
                app.sldr_full.Limits = [full_idx(1), full_idx(2)-obj.dur_ov+2];
            else
                app.sldr_full.Limits = [full_idx(1), full_idx(2)-obj.dur_ov+1];
                app.sldr_full.Enable = 'on';
            end
            app.sldr_overview.Limits(1) = 1;
            app.sldr_overview.Limits(2) = obj.data_len;
            app.sldr_overview.Limits(1) = obj.start_idx_ov;
            app.sldr_overview.Limits(2) = obj.stop_idx_ov-obj.dur_det+1;
            
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

            notify (obj, 'new_frame_full')

        end


        function obj = set_ov_dur(obj,dur,app)
            obj.dur_ov = int64(dur*1/mean(diff(app.data(:,2))));
            
            tmp = [obj.start_idx_full, obj.stop_idx_full-obj.dur_ov+1];
            if tmp(1) == tmp(2)
                app.sldr_full.Enable = 'off';
                app.sldr_full.Limits(1) = 1;
                app.sldr_full.Limits(2) = obj.data_len;
                app.sldr_full.Limits(1) = obj.start_idx_full;
                app.sldr_full.Limits(2) = obj.stop_idx_full-obj.dur_ov+2;
            else
                app.sldr_full.Limits(1) = 1;
                app.sldr_full.Limits(2) = obj.data_len;
                app.sldr_full.Limits(1) = obj.start_idx_full;
                app.sldr_full.Limits(2) = obj.stop_idx_full-obj.dur_ov+1;
                app.sldr_full.Enable = 'on';
            end
            
            app.sldr_overview.Limits(1) = 1;
            app.sldr_overview.Limits(2) = obj.data_len;
            app.sldr_overview.Limits(1) = obj.start_idx_full;
            app.sldr_overview.Limits(2) = obj.stop_idx_full-obj.dur_ov+1;
            obj.stop_idx_ov = obj.start_idx_ov+obj.dur_ov;
            if obj.start_idx_ov < 1
                obj.start_idx_ov = int64(1);
                obj.stop_idx_ov = int64(obj.start_idx_ov +obj.dur_ov);
            end
            if obj.stop_idx_ov > obj.stop_idx_full
                obj.stop_idx_ov = int64(obj.stop_idx_full);
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
            app.sldr_overview.Limits(1) = 1;
            app.sldr_overview.Limits(2) = obj.data_len;
            app.sldr_overview.Limits(1) = obj.start_idx_ov;
            app.sldr_overview.Limits(2) = obj.stop_idx_ov-obj.dur_det+1;
            
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

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

            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);
            
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

            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

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

            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

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
            
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

            notify(obj, 'new_frame')
        end

        function obj = get_frame_ov(obj, idx,app)
            obj.start_idx_ov = int64(idx);
            obj.stop_idx_ov = idx +obj.dur_ov;
            if obj.start_idx_ov < 1
                obj.start_idx_ov = int64(1);
                obj.stop_idx_ov = int64(obj.start_idx_ov +obj.dur_ov);
            end
            if obj.stop_idx_ov > obj.stop_idx_full
                obj.stop_idx_ov = int64(obj.stop_idx_full);
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
            
            spk_idx = find((app.spike_res.spike_idx>= obj.start_idx_ov) & (app.spike_res.spike_idx<= obj.stop_idx_ov));
            displ_spks =  [app.spike_res.spike_idx(spk_idx)',app.spike_res.use_spikes(spk_idx,1),  app.spike_res.cluster(spk_idx)'];
            obj.ov_idx = displ_spks(logical(displ_spks(:,2)),1);

            notify(obj, 'new_frame_ov')
        end
            
    end
    
end

