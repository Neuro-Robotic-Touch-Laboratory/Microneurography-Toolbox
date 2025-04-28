function create_mat_file(filename,varargin)
%create_mat_file function to create .mat files readable by the MNG toolbox
%
%   create_mat_file(filename,varargin) create a mat file named filename
%   defined by varargin
%   create_mat_file(__, Name, Values)
%   
%   Name-Values Arguments
%
%   'channel' - define channel with values:
%       data : array with datapoints
%       timestamps / sampling frequency, start : either array with
%           timestamps same size as data or array with length 2 containing
%           sampling frequency and start timestamps 
%       name : char array containing channel name 
%       units : char array containing channel units
%
%   'comments' - define comments
%       comment strings : cell array containing the comments
%       timestamps : array containing the timestamps 


channel_version = 1;
comment_version = 1;
data_version = 1;
file_version = 1;
record_version = 1;

temp_data = struct();
data_idx=1;
temp_comm = struct;

while ~isempty(varargin)
    switch varargin{1}
        case 'channel'
            data = varargin{2};
            ts_fs = varargin{3};
            name = varargin{4};
            unit = varargin{5};

            if size(data,1) > size(data,2)
                temp_data(data_idx).data = data;
            else
                temp_data(data_idx).data = data';
            end
            temp_data(data_idx).n_samples = length(data);

            if length(ts_fs) == 2 
                temp_data(data_idx).ts = ((0:length(data)-1)/ts_fs(1)+ts_fs(2))';
            else
                if size(ts_fs,1) > size(ts_fs,2)
                    temp_data(data_idx).ts = ts_fs;
                else
                    temp_data(data_idx).ts = ts_fs';
                end
            end
            temp_data(data_idx).dt = round(mean(diff(temp_data(end).ts)),6);
            temp_data(data_idx).name = name;

            temp_data(data_idx).units = unit;
            
            data_idx = data_idx+1;
            clear data ts_fs name unit
            varargin(1:5) = [];

        case 'comments'
            comms = varargin{2};
            ts = varargin{3};
            for i = 1: length(comms)
                temp_comm(i).com_str = comms{i};
                temp_comm(i).ts = ts(i);
            end
            clear comms ts
            varargin(1:3) = [];

        otherwise
            varargin(1) = [];
    end
        
end 

clear varargin data_idx

[~,idx] = sort(vertcat(temp_data.dt));
file_meta = struct();
file_meta.n_records = 1;
file_meta.n_channels = length(idx);
channel_meta = struct();
for i =1 : length(idx)
    channel_meta(i).id = i;
    channel_meta(i).name = temp_data(idx(i)).name;
    channel_meta(i).units = {temp_data(idx(i)).units};
    channel_meta(i).n_samples = temp_data(idx(i)).n_samples;
    channel_meta(i).dt = temp_data(idx(i)).dt;
    eval(['data__chan_' num2str(i) '_rec_1 = temp_data(idx(i)).data;'])
end
record_meta = struct();
record_meta.n_ticks = length(temp_data(idx(1)).data);
record_meta.tick_dt = temp_data(idx(1)).dt;
record_meta.record_start = temp_data(idx(1)).ts(1);
record_meta.data_start = temp_data(idx(1)).ts(1);
record_meta.trigger_minus_rec_start_samples = 0;

comments = struct();

for i = 1 :length(temp_comm)
    comments(i).str = temp_comm(i).com_str;
    comments(i).id = i;
    comments(i).tick_position = find(temp_data(idx(1)).ts >= temp_comm(i).ts,1,'first');
    comments(i).channel = -1;
    comments(i).record = 1;
    comments(i).tick_dt = temp_data(idx(1)).dt;
end
clear idx temp_data temp_comm i

save(filename)
disp([filename '.mat was successfully created'])
disp(['Start time: ' num2str(record_meta.data_start) ' s, data length: ' num2str(record_meta.n_ticks*record_meta.tick_dt) ' s'])
disp(['the file contains ' num2str(length(channel_meta)) ' channels:'])

for i = 1 : length(channel_meta)
    disp(['Channel ' num2str(i) ' name: ' channel_meta(i).name ', units: ' channel_meta(i).units{1} ', length: ' num2str(channel_meta(i).n_samples) ' samples, sampling frequency: ' num2str(1/channel_meta(i).dt) ' Hz'])
end

disp(['the file contains ' num2str(length(comments)) ' comments:'])
for i = 1 : length(comments)
    disp(['Comment ' num2str(i) ': ' comments(i).str ', at: ' num2str((comments(i).tick_position-1)*comments(i).tick_dt +record_meta.data_start) ' s'])
end

end