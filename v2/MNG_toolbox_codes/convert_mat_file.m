function file_path_out = convert_mat_file(file_path_in)
dot_idx = find (char(file_path_in)=='.',1,'last');
file_path_out = [file_path_in(1:dot_idx-1) '_tmp.mat' ];
file_in = load(file_path_in);

file_meta.n_records = size(file_in.datastart,2);
file_meta.n_channels = size(file_in.datastart,1);
file_version = 1;
for i = 1: file_meta.n_records
    record_meta(i).tick_dt = 1/ file_in.tickrate(i);
    record_meta(i).record_start = file_in.blocktimes(i);
    record_meta(i).data_start = file_in.blocktimes(i);
    record_meta(i).trigger_minus_rec_start_samples = 0;
end
record_version = 1;
comments = struct('str', [], 'id', [], 'tick_position', [], 'channel', [], 'record', [], 'tick_dt',[]);
comments(1) = [];
for i = 1: size(file_in.com,1)
    comments(i).str = strtrim(file_in.comtext(file_in.com(i,5),:));
    comments(i).id = file_in.com(i,5);
    comments(i).tick_position = file_in.com(i,3);
    comments(i).channel = file_in.com(i,1);
    comments(i).record = file_in.com(i,2);
    comments(i).tick_dt = record_meta(file_in.com(i,2)).tick_dt;
end

channel_meta.id = [];
channel_meta.name =  [];
channel_meta.units = {[]};
channel_meta.n_samples = [];
channel_meta.dt = [];

channel_version = 1;
comment_version = 1;

[~,min_dt_idx] = max (file_in.samplerate(:,1));

for i = 1 : file_meta.n_records
    for j = 1 : file_meta.n_channels
        name = ['data__chan_' num2str(j), '_rec_' num2str(i)];
        data = file_in.data(file_in.datastart(j,i):file_in.dataend(j,i))';
        create_var(name,data)
        if i ==1
            channel_meta(j).id = j;
            channel_meta(j).name =  strtrim(file_in.titles(j,:));
        end
        channel_meta(j).n_samples(i) = length(data);
        channel_meta(j).dt(i) = 1/ file_in.samplerate(j,i);
        if file_in.unittextmap(i) > 0 
            channel_meta(j).units{i} = strtrim(file_in.unittext(file_in.unittextmap(i),:));
        else 
            channel_meta(j).units{i} = '';
        end
        
    end
    record_meta(i).tick_dt = 1/ file_in.tickrate(1);
    record_meta(i).record_start = file_in.blocktimes;
    record_meta(i).data_start = file_in.blocktimes;
    record_meta(i).trigger_minus_rec_start_samples = 0;
    record_meta(i).n_ticks = channel_meta(min_dt_idx).n_samples(i)*(record_meta(i).tick_dt/channel_meta(min_dt_idx).dt(i));
    
end


clear dot_idx file_path_in i j name data file_in min_dt_idx
save(file_path_out,'-v7.3')
end

function create_var(name,data)
    assignin('caller',name,data)
end



