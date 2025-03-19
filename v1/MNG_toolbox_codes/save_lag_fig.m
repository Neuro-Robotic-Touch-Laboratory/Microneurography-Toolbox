function save_lag_fig(app)


print_cell = {};
print_cell{1,1} = 'Interval';
print_cell{1,2} = 'signal1';
print_cell{1,3} = 'movmean';
print_cell{1,4} = 'signal2';
print_cell{1,5} = 'movmean';
print_cell{1,6} = 'Lag';
print_cell{1,7} = 'Kendalls tau';
print_cell{1,8} = 'p';
print_cell{1,9} = 'spearmans rho';
print_cell{1,10} = 'p';

for i=1:length(app.lag)
    app.popup_lag_bgsig.Value = app.popup_lag_bgsig.Items{app.lag(i).bgsig};
    app.popup_lag_int.Value = app.popup_lag_int.Items{app.lag(i).int};
    app.popup_lag_sig1.Value = app.popup_lag_sig1.Items{app.lag(i).sig1};
    app.popup_lag_sig2.Value = app.popup_lag_sig2.Items{app.lag(i).sig2};
    if app.lag(i).mm1 == 0
        app.chkbx_lag_ma_sig1.Value = false;
    else
        app.chkbx_lag_ma_sig1.Value = true;
        app.edt_lag_ma_sig1.Value = app.lag(i).mm1;
    end

    if app.lag(i).mm2 == 0
        app.chkbx_lag_ma_sig2.Value = false;
    else
        app.chkbx_lag_ma_sig2.Value = true;
        app.edt_lag_ma_sig2.Value = app.lag(i).mm2;
    end
    update_lag_axis(app, 'sig')
    update_lag_axis(app, 'calc')
    update_lag_axis(app, 'select',app.lag(i).lag)
    app.sldr_lag_start.Value = app.lag(i).lims(1);
    app.sldr_lag_end.Value = app.lag(i).lims(2);
    [rk,pk, rs,ps] = plot_stuff(app, app.lag(i));
    
    print_cell{end+1,1} = app.popup_lag_int.Value;
    print_cell{end,2} = app.popup_lag_sig1.Value;
    print_cell{end,3} = app.lag(i).mm1;
    print_cell{end,4} = app.popup_lag_sig2.Value;
    print_cell{end,5} = app.lag(i).mm2;
    print_cell{end,6} = app.lag(i).lag;
    print_cell{end,7} = rk;
    print_cell{end,8} = pk;
    print_cell{end,9} = rs;
    print_cell{end,10} = ps;
    
    disp(['i: ' num2str(i) ' / ' num2str(length(ints)) ' , j: ' num2str(j) ' / ' num2str(length(sig1)) ' , k: ' num2str(k)  '/ ' num2str(length(sig2))])
    
end

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';
writecell(print_cell,[app.settings.output_dir '\' file '_INT_' num2str(round(app.settings.interval(1,1),1)) '-' num2str(round(app.settings.interval(1,2),1)) 'spike_results.xls'])

end

function [rk,pk, rs,ps] = plot_stuff(app, lag_set)
%%
tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
h = figure;
h.Position = [20,100,1400,840];

subplot(3,5,1:3)

mn1 = min(app.settings.lagsig1);
mx1 = max(app.settings.lagsig1);
if ~isnan(app.settings.lagbgsig(1,1))
    plot (app.settings.lagbgsig(:,2),app.settings.lagbgsig(:,1)*diff([mn1,mx1])+mn1, Color=[0.7 0.7 0.7])
else
    plot (nan,nan, Color=[0.7 0.7 0.7])
end
hold on
plot (app.settings.lagts,app.settings.lagsig1, Color=[0 0.4470 0.7410])
plot (app.settings.lagts(1:app.sldr_lag_start.Value),app.settings.lagsig1(1:app.sldr_lag_start.Value), Color=[0.8,0.8,0.8]);
plot (app.settings.lagts(app.sldr_lag_end.Value:length(app.settings.lagts)),app.settings.lagsig1(app.sldr_lag_end.Value:length(app.settings.lagts)), Color=[0.8,0.8,0.8]);
hold off
xlim ([app.settings.lagts(1),app.settings.lagts(end)])



if app.chkbx_lag_ma_sig1.Value
    tmp_str = [num2str(app.edt_lag_ma_sig1.Value) ' s movmean'];
else
    tmp_str = '';
end

title ([app.popup_lag_sig1.Value ' ' tmp_str])

subplot(3,5,6:8)

mn2 = min(app.settings.lagsig2);
mx2 = max(app.settings.lagsig2);
if ~isnan(app.settings.lagbgsig(1,1))
    plot (app.settings.lagbgsig(:,2),app.settings.lagbgsig(:,1)*diff([mn2,mx2])+mn2, Color=[0.7 0.7 0.7])
else
    plot (nan,nan, Color=[0.7 0.7 0.7])
end
hold on
plot (app.settings.lagts,app.settings.lagsig2, Color=[0.8500 0.3250 0.0980])
plot (app.settings.lagts(1:app.sldr_lag_start.Value),app.settings.lagsig2(1:app.sldr_lag_start.Value), Color=[0.8,0.8,0.8]);
plot (app.settings.lagts(app.sldr_lag_end.Value:length(app.settings.lagts)),app.settings.lagsig2(app.sldr_lag_end.Value:length(app.settings.lagts)), Color=[0.8,0.8,0.8]);
hold off
xlim ([app.settings.lagts(1),app.settings.lagts(end)])

if app.chkbx_lag_ma_sig2.Value
    tmp_str = [num2str(app.edt_lag_ma_sig2.Value) ' s movmean'];
else
    tmp_str = '';
end

title ([app.popup_lag_sig2.Value ' ' tmp_str])


subplot(3,5,[4,5,9,10])

plot (app.settings.lagxcts,app.settings.lagxc)
xlim ([app.settings.lagxcts(1),app.settings.lagxcts(end)])
line([lag_set.lag,lag_set.lag],ylim,'Color', 'r')
title (['Lag: ' num2str(lag_set.lag) 's'])
subplot(3,5,11:15)

tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];

sig1 = app.settings.lagsig1(tmp(1):tmp(2))-mean(app.settings.lagsig1(tmp(1):tmp(2)));
sig1 = sig1/max(abs(sig1));

sig2 = app.settings.lagsig2(tmp(1):tmp(2))-mean(app.settings.lagsig2(tmp(1):tmp(2)));
sig2 = sig2/max(abs(sig2));
idx =  lag_set.lag/mean(diff(app.settings.lagts)); 
if lag_set.lag>=0
    sig2 = [nan(idx,1);sig2];
else
    sig1 = [nan(abs(idx),1);sig1];
end
ts1 = (1:length(sig1))*mean(diff(app.settings.lagts));
ts2 = (1:length(sig2))*mean(diff(app.settings.lagts));
if app.settings.inverse_lagsig
    sig2 = sig2*(-1);
end

plot(ts1,sig1, LineWidth=1.5)
hold ('on') 
plot(ts2,sig2, LineWidth=1.5)
hold ('off')
xlim ([ts1(1),max(ts1(end),ts2(end))])

sig1 = app.settings.lagsig1(tmp(1):tmp(2))-mean(app.settings.lagsig1(tmp(1):tmp(2)));
sig1 = sig1/max(abs(sig1));

sig2 = app.settings.lagsig2(tmp(1):tmp(2))-mean(app.settings.lagsig2(tmp(1):tmp(2)));
sig2 = sig2/max(abs(sig2)); 
if idx>=0
    sig1(1:idx) = [];
    sig2(end-idx+1:end) = [];
else
    sig2(1:idx) = [];
    sig1(end-idx+1:end) = [];
end
rk=999;pk =999;
%     [rk,pk ]  = corr(sig1,sig2,'Type','Kendall');
[rs,ps ]  = corr(sig1,sig2,'Type','Spearman');

title(['Spearmans r: ' num2str(rs) ', p: ' num2str(rs)])


bs = find (app.settings.file_path == '\',1,"last");
dt = find (app.settings.file_path == '.',1,"last");

if app.chkbx_lag_ma_sig1.Value
    mm1  = [ '_ma-' num2str(app.edt_lag_ma_sig1.Value)];
else
    mm1  = '';
end

if app.chkbx_lag_ma_sig2.Value
    mm2  = [ '_ma-' num2str(app.edt_lag_ma_sig2.Value)];
else
    mm2  = '';
end

fn = [app.settings.output_dir '\' app.settings.file_path(bs+1:dt-1) '-' app.popup_lag_int.Value '-' app.popup_lag_sig1.Value mm1 '-' app.popup_lag_sig2.Value mm2];
savefig(h,fn)
close(h)
end

