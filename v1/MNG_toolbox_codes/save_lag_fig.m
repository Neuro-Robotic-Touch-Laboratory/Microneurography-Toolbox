function save_lag_fig(app)
[ints,~] = listdlg('ListString',app.popup_lag_int.Items, 'PromptString','Select a intervals');
[sig1, ~]  = listdlg('ListString',app.popup_lag_sig1.Items, 'PromptString','Select a signal 1');
[sig2, ~]  = listdlg('ListString',app.popup_lag_sig2.Items, 'PromptString','Select a signal 2');

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

for i=1:length(ints)
    app.popup_lag_int.Value = app.popup_lag_int.Items{ints(i)};
    for j = 1:length(sig1)
        app.popup_lag_sig1.Value = app.popup_lag_sig1.Items(sig1(j));
        for k =1:length(sig2)
            app.popup_lag_sig2.Value = app.popup_lag_sig2.Items(sig2(k));
%             app.chkbx_lag_ma_sig1.Value = false;
%             app.chkbx_lag_ma_sig2.Value = false;
% %             app.chkbx_invert_sig1.Value = false;
% %             app.chkbx_invert_sig2.Value = false;
%             update_lag_axis(app, 'sig')
%             update_lag_axis(app, 'calc')
%             [lag,rk,pk, rs,ps] = plot_stuff(app);
%             
%             print_cell{end+1,1} = app.popup_lag_int.Value;
%             print_cell{end,2} = app.popup_lag_sig1.Value;
%             print_cell{end,3} = '0';
%             print_cell{end,4} = app.popup_lag_sig2.Value;
%             print_cell{end,5} = 0;
%             print_cell{end,6} = lag;
%             print_cell{end,7} = rk;
%             print_cell{end,8} = pk;
%             print_cell{end,9} = rs;
%             print_cell{end,10} = ps;

%             app.chkbx_lag_ma_sig1.Value = true;
%             app.chkbx_lag_ma_sig2.Value = true;
%             app.edt_lag_ma_sig1.Value = 1;
%             app.edt_lag_ma_sig2.Value = 1;
% 
%             update_lag_axis(app, 'sig')
%             update_lag_axis(app, 'calc')
%             [lag,rk,pk, rs,ps] = plot_stuff(app);
% 
%             print_cell{end+1,1} = app.popup_lag_int.Value;
%             print_cell{end,2} = app.popup_lag_sig1.Value;
%             print_cell{end,3} = 1;
%             print_cell{end,4} = app.popup_lag_sig2.Value;
%             print_cell{end,5} = 1;
%             print_cell{end,6} = lag;
%             print_cell{end,7} = rk;
%             print_cell{end,8} = pk;
%             print_cell{end,9} = rs;
%             print_cell{end,10} = ps;

            app.chkbx_lag_ma_sig1.Value = true;
            app.chkbx_lag_ma_sig2.Value = true;
            app.edt_lag_ma_sig1.Value = 3;
            app.edt_lag_ma_sig2.Value = 3;

            update_lag_axis(app, 'sig')
            update_lag_axis(app, 'calc')
            [lag,rk,pk, rs,ps] = plot_stuff(app);

            print_cell{end+1,1} = app.popup_lag_int.Value;
            print_cell{end,2} = app.popup_lag_sig1.Value;
            print_cell{end,3} = 3;
            print_cell{end,4} = app.popup_lag_sig2.Value;
            print_cell{end,5} = 3;
            print_cell{end,6} = lag;
            print_cell{end,7} = rk;
            print_cell{end,8} = pk;
            print_cell{end,9} = rs;
            print_cell{end,10} = ps;

            
            app.chkbx_lag_ma_sig1.Value = true;
            app.chkbx_lag_ma_sig2.Value = true;
            app.edt_lag_ma_sig1.Value = 5;
            app.edt_lag_ma_sig2.Value = 5;
%             app.chkbx_invert_sig1.Value = false;
%             app.chkbx_invert_sig2.Value = true;
            update_lag_axis(app, 'sig')
            update_lag_axis(app, 'calc')
            [lag,rk,pk, rs,ps] = plot_stuff(app);
            
            print_cell{end+1,1} = app.popup_lag_int.Value;
            print_cell{end,2} = app.popup_lag_sig1.Value;
            print_cell{end,3} = 5;
            print_cell{end,4} = app.popup_lag_sig2.Value;
            print_cell{end,5} = 5;
            print_cell{end,6} = lag;
            print_cell{end,7} = rk;
            print_cell{end,8} = pk;
            print_cell{end,9} = rs;
            print_cell{end,10} = ps;

            disp(['i: ' num2str(i) ' / ' num2str(length(ints)) ' , j: ' num2str(j) ' / ' num2str(length(sig1)) ' , k: ' num2str(k)  '/ ' num2str(length(sig2))])
        end
    end
end

writecell(print_cell,[app.settings.output_dir '\sig_lags.xls' ])

end

function [lag,rk,pk, rs,ps] = plot_stuff(app)
tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
h = figure;
h.Position = [20,100,1120,840];

subplot(2,2,2)
tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
plot(app.settings.lagxcts, app.settings.lagxc)
[~,idx] = max(abs(app.settings.lagxc));
line([app.settings.lagxcts(idx), app.settings.lagxcts(idx)], ylim, 'Color', 'r')
line([app.settings.lagxcts(1) app.settings.lagxcts(end)], [0,0], 'Color', 'k','LineStyle',':')

lag = (idx-length(app.settings.lagsig1(tmp(1):tmp(2))))*mean(diff(app.settings.lagts));

idx = idx-length(app.settings.lagsig1(tmp(1):tmp(2)));
xlim ([app.settings.lagxcts(1), app.settings.lagxcts(end)])

ax = subplot(2,2,1);
ax.XTick = [0,1];
ax.YTick = [];
ylim ([0,12])
xlim ([0,2 ])
text(0.1,10,['Intervall: ' char(app.popup_lag_int.Value)])
if app.chkbx_lag_ma_sig1.Value
tmp_str =[ ' movmean: ' num2str(app.edt_lag_ma_sig1.Value) ' s'];
else
tmp_str = '';
end
text (0.1,8.5,['Signal 1: ' (app.popup_lag_sig1.Value) tmp_str])

if app.chkbx_lag_ma_sig2.Value
tmp_str =[ ' movmean: ' num2str(app.edt_lag_ma_sig2.Value) ' s'];
else
tmp_str = '';
end
text(0.1,7,['Signal 1: ' (app.popup_lag_sig2.Value) tmp_str])
text(0.1,5.5,['Lag: ' num2str(lag) ' s'])
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
% [rk,pk ]  = corr(sig1,sig2,'Type','Kendall');
[rs,ps ]  = corr(sig1,sig2,'Type','Spearman');

if app.settings.inverse_lagsig
    text(0.1,4, 'signal 2 inversed')
else
    text(0.1,4, 'signal 2 not inversed')
end
text(0.1,2.5,['Spearmans rho: ' num2str(rs,3) ', p: ' num2str(ps,5)])

subplot(2,2,3:4)

tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];

sig1 = app.settings.lagsig1(tmp(1):tmp(2))-mean(app.settings.lagsig1(tmp(1):tmp(2)));
sig1 = sig1/max(abs(sig1));

sig2 = app.settings.lagsig2(tmp(1):tmp(2))-mean(app.settings.lagsig2(tmp(1):tmp(2)));
sig2 = sig2/max(abs(sig2));
if idx>=0
    sig2 = [nan(idx,1);sig2];
else
    sig1 = [nan(abs(idx),1);sig1];
end
if app.settings.inverse_lagsig
    sig2 = sig2*(-1);
end

plot(sig1, LineWidth=1.5)
hold ('on') 
plot(sig2, LineWidth=1.5)
hold ('off')
xlim ([1,max(length(sig1),length(sig2))])

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