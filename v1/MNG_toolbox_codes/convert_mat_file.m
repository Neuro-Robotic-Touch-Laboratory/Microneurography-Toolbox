function file_path_out = convert_mat_file(file_path_in)
dot_idx = find (char(file_path_in)=='.',1,'last');
file_path_out = [file_path_in(1:dot_idx-1) '_tmp.mat' ];
file_in = load(file_path_in);

file_meta.n_records = size(file_in.datastart,2);
file_meta.n_channels = size(file_in.datastart,1);
file_version = 1;

record_meta.tick_dt = 1/ file_in.tickrate(1);
record_meta.record_start = file_in.blocktimes;
record_meta.data_start = file_in.blocktimes;
record_meta.trigger_minus_rec_start_samples = nan;
record_version = 1;

for i = 1: length(file_in.com)
    comments(i).str = strtrim(file_in.comtext(file_in.com(i,5),:));
    comments(i).id = i;
    comments(i).tick_position = file_in.com(i,3);
    comments(i).channel = file_in.com(i,1);
    comments(i).record = file_in.com(i,2);
    comments(i).tick_dt = record_meta.tick_dt;
end

for i = 1 : file_meta.n_channels
    channel_meta(i).id = i;
    channel_meta(i).name =  strtrim(file_in.titles(i,:));
    if file_in.unittextmap(i) > 0 
        channel_meta(i).units = {strtrim(file_in.unittext(file_in.unittextmap(i),:))};
    else 
        channel_meta(i).units = {''};
    end
    channel_meta(i).n_samples = 0;
    channel_meta(i).dt = 1/ file_in.samplerate(i,1);
end
channel_version = 1;
comment_version = 1;


for i = 1 : file_meta.n_records
    for j = 1 : file_meta.n_channels
        name = ['data__chan_' num2str(j), '_rec_' num2str(i)];
        data = file_in.data(file_in.datastart(j,i):file_in.dataend(j,i))';
        create_var(name,data)
        channel_meta(j).n_samples = channel_meta(j).n_samples +length(data);
    end
end
[~,idx] = max(file_in.samplerate); 
record_meta.n_ticks = channel_meta(idx).n_samples; 
clear dot_idx file_path_in i j name data
save(file_path_out,'-v7.3')
end

function create_var(name,data)
    assignin('caller',name,data)
end
