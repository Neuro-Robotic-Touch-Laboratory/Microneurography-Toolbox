function reset_results(app,res_data,res_derived)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if res_data
    app.data = [];
    app.comments = [];
    app.settings = struct('interval', [nan, nan; nan, nan],...
                          'channel_sheet', [nan, nan],...
                          'click_select', [false, false],...
                          'file_path','',...
                          'done', struct('hb', false,...
                                         'intervals', false,...
                                         'burst', false, ...
                                         'resp',false),...
                          'burst_mod_mode', [0,nan],...
                          'burst_rem_int',[], ...
                          'channel_idx', struct('msna',nan, 'ecg', nan, 'bldp', nan, 'resp', nan),...
                          'output_dir', 'none', ...
                          'figs', [],...
                          'rep_dir', 'none',...
                          'brst_files',[],...
                          'int_templ',[],...
                          'collect_path','none',...
                          'tick_freq', nan...
                          );
    app.hb_res = [];
    app.burst_ints = [];
    app.burst_res = [];
    app.spike_res = [];
    app.spec_res = [];
    app.entropy_res = [];
    app.stepwise_res = [];
    app.resp_res = []; 
    app.bp_res = [];
    app.corre_res = [];
    app.lbl_output.Text = 'none';
    app.comb_stat = [];
end

if res_derived
    rem_idx = [];
    for i =1:length(app.data)
        if app.data(i).derived
            rem_idx = [rem_idx,i];
        end
    end
    if ~isempty(rem_idx)
        app.data(rem_idx) = [];
    end
    app.hb_res = [];
    app.burst_ints = [];
    app.burst_res = [];
    app.spike_res = [];
    app.spec_res = [];
    app.entropy_res = [];
    app.stepwise_res = [];
    app.resp_res = []; 
    app.bp_res = [];
    app.corre_res = [];
    app.settings.done.hb = false;
    app.settings.done.intervals = false;
    app.settings.done.burst = false;
    app.settings.done.resp = false;
%     app.settings.interval(1,:) = [1, app.settings.interval(2,2)];
end
app.btn_hb.BackgroundColor = [.96, .96, .96];
app.btn_sel_int.BackgroundColor = [.96, .96, .96];
app.btn_resp.BackgroundColor = [.96, .96, .96];
app.btn_burst_analysis.BackgroundColor = [.96, .96, .96];
app.popup_int_spect.Items = {'full'};
update_intervall_popup(app)
update_signal_popup(app)
clear_axis(app)
end

