function burst_res = batch_bursts(data,ts_fs,varargin)
app = struct();
app.edt_window.Value = 0.18;
app.chkbx_detrent_rms.Value = false;
app.edt_threshold.Value = 0.11;
app.edt_threshold_abs.Value = 0.11;
app.chkbx_manual.Value = false;
app.chkbx_th_on_silent.Value = false;
app.edt_outlier_threshold.Value = 100;
while ~isempty(varargin)
    switch varargin{1}
        case 'window'
        	app.edt_window.Value = varargin{2};
            varargin(1:2) = [];
        case 'threshold'
        	app.edt_threshold.Value = varargin{2};
            varargin(1:2) = [];
        case 'detrent'
            app.chkbx_detrent_rms.Value = varargin{2};
            varargin(1:2) = [];
        case 'manual'
            app.edt_threshold_abs.Value = varargin{2};
            app.chkbx_manual.Value = true;
            varargin(1:2) = [];
        case 'onsilent'
            app.chkbx_th_on_silent.Value = varargin{2};
            varargin(1:2) = [];
        case 'outlier'
            app.edt_outlier_threshold.Value = varargin{2};
            varargin(1:2) = [];
        otherwise
            varargin(1) = [];
    end
end
if length(ts_fs) == length(data)
    ts = ts_fs;
    fs = 1/mean(diff(ts));
else
    fs = ts_fs;
    ts = (1:length(data))/fs;
end

app.settings.interval = [ts(1), ts(end)];
app.settings.channel_idx.msna = 1;
app.data(1).data = data;
app.data(1).ts(1)= 1/fs;
app.data(1).unit = '[n/a]';



burst_res = burst_analysis(app);