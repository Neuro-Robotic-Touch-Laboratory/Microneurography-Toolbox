 function save_session (app)
datafile = app.settings.file_path;
savestruct = struct();
tmp = find(datafile == '\',1,"last");
[file, path] = uiputfile ('*.ssf' , 'select file / path to save the session', [datafile(tmp:end) '_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_save.ssf'] );
savestruct.settings.datafile = datafile;
channels = {};

for i = 1 : length(app.data)
    if ~app.data(i).derived & ~strcmp(app.data(i).name, 'RMS MSNA')
        channels{end+1} = app.data(i).name;
    end
end
savestruct.settings.channels = channels;
savestruct.settings.channel_idx = app.settings.channel_idx;
savestruct.settings.interval =  app.settings.interval;
savestruct.settings.res_folder = app.settings.output_dir;
savestruct.settings.int_templ = app.settings.int_templ;
if app.settings.done.hb
    savestruct.hb.idx = app.hb_res.idx;
    savestruct.hb.ts = app.hb_res.t_events;
end

if app.settings.done.resp
    savestruct.resp.idx = app.resp_res.idx;
    savestruct.resp.ts = app.resp_res.ts;
end

if app.settings.done.intervals
    savestruct.intervals.idx = app.burst_ints;
end

if app.settings.done.burst
    savestruct.burst.detrend =  app.chkbx_detrent_rms.Value;

%     savestruct.burst.outlier_thresh =  app.edt_outlier_threshold.Value;
%     savestruct.burst.window =  app.edt_window.Value;
%     savestruct.burst.thresh =  app.edt_threshold.Value;
    savestruct.burst.window = app.burst_res.window;
    savestruct.burst.delay = app.burst_res.delay;
    savestruct.burst.use_burst =  app.burst_res.use_burst;
    savestruct.burst.rem_int =  app.settings.burst_rem_int;  
end

if ~isempty(app.spike_res)
    savestruct.spike.spike_idx = app.spike_res.spike_idx;
    savestruct.spike.spike_ts = app.spike_res.spike_ts;
    savestruct.spike.cluster = app.spike_res.cluster;
    savestruct.spike.extremes = app.spike_res.extremes;
    savestruct.spike.use_spikes = app.spike_res.use_spikes;
    savestruct.spike.analysis = app.spike_res.analysis  ;
    savestruct.spike.par = app.settings.par;
    savestruct.spike.spike = app.spike_res.spike;
end

save( [path  file],'savestruct','-mat')


