function update_lag_axis(app, action,varargin)

load = false;
plt = false;
update_sig = false; 
update_sldr = false;
update_ov = false;
calc_corr = false;
app.settings.inverse_lagsig = false;
switch action
    case 'init'

        load = true;
        plt = true;
    case 'sig'
        
        load = true;
        update_sig = true;
        update_sldr = true;

    case 'slider'
        update_sldr = true;
        cla(app.ax_lag_xcorr)
        cla(app.ax_lag_overlay)
    case 'calc'
        tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
        sig1 = app.settings.lagsig1(tmp(1):tmp(2))-mean(app.settings.lagsig1(tmp(1):tmp(2)));
        sig1 = sig1/max(abs(sig1));
        
        sig2 = app.settings.lagsig2(tmp(1):tmp(2))-mean(app.settings.lagsig2(tmp(1):tmp(2)));
        sig2 = sig2/max(abs(sig2));

        [app.settings.lagxcts, app.settings.lagxc, lag] = calc_lac(sig1,sig2,app);
        
        plot(app.ax_lag_xcorr,app.settings.lagxcts, app.settings.lagxc,'HitTest','off')

        line(app.ax_lag_xcorr,[app.settings.lagxcts(1) app.settings.lagxcts(end)], [0,0], 'Color', 'k','LineStyle',':','HitTest','off')
        line(app.ax_lag_xcorr,[lag lag], ylim(app.ax_lag_xcorr), 'Color', 'r','HitTest','off')
%         lag = (idx-length(app.settings.lagsig1(tmp(1):tmp(2))))*mean(diff(app.settings.lagts));
        app.lbl_lag_lag.Text = ['Lag: ' num2str(lag) ' s'];
%         idx = idx-length(app.settings.lagsig1(tmp(1):tmp(2)));
        xlim (app.ax_lag_xcorr, [app.settings.lagxcts(1), app.settings.lagxcts(end)])
        update_ov = true;
        calc_corr = true;
        app.settings.lag = lag;

    case 'select'
%         loc = round(varargin{1,1});
        loc = find (app.settings.lagxcts >= varargin{1,1},1);
        tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
        if app.chkbx_lock.Value
            idx= [];
            rng = round(length(app.settings.lagxcts)/40);
            len = length(app.settings.lagxcts);
            while isempty(idx)    
                if loc-rng<1 | loc+rng>len
                    return
                end

                if app.chkbx_lag_wavelet.Value
                    tmp_int = app.settings.lagxc(loc-rng: loc+rng);
                else
                    tmp_int = abs(app.settings.lagxc(loc-rng: loc+rng));
                end
                
                [~,idx] = findpeaks(tmp_int, 'NPeaks',1);
    
                idx = round(idx-1-rng+loc);
                if isempty(idx)
                    if tmp_int(1) <tmp_int(end)
                        step = round(1.5*rng);
                    else
                        step = round(1.5*rng)*(-1);
                    end
                    loc = loc+step;
                end
            end

        else
            rng = round(0.5/mean(diff(app.settings.lagts)));
            [~,idx] = max(abs(app.settings.lagxc(loc-rng: loc+rng)));
            
            idx = round(idx-1-rng+loc);
        end

        if app.settings.lagxc(idx) <0 
            app.settings.inverse_lagsig = true;
        end

        plot(app.ax_lag_xcorr, app.settings.lagxcts, app.settings.lagxc,'HitTest','off')
        line(app.ax_lag_xcorr,[app.settings.lagxcts(idx), app.settings.lagxcts(idx)], ylim(app.ax_lag_xcorr), 'Color', 'r','HitTest','off')
        line(app.ax_lag_xcorr,[app.settings.lagxcts(1) app.settings.lagxcts(end)], [0,0], 'Color', 'k','LineStyle',':','HitTest','off')
        lag = app.settings.lagxcts(idx);
        app.lbl_lag_lag.Text = ['Lag: ' num2str(lag) ' s'];
        idx = idx-length(app.settings.lagsig1(tmp(1):tmp(2)));
        xlim (app.ax_lag_xcorr, [app.settings.lagxcts(1),app.settings.lagxcts(end)])
        app.settings.lag = lag;
        update_ov = true;
        calc_corr = true;
end

if load 
        tmp = find(strcmp(app.popup_lag_sig1.Value, app.popup_lag_sig1.Items));
        [tmp_data1,tmp_ts1,tmp_name1, tmp_unit1] = current_signal(app, tmp);

        tmp_data1(:,2)= tmp_ts1(1):tmp_ts1(1):tmp_ts1(2);

        tmp = find(strcmp(app.popup_lag_sig2.Value, app.popup_lag_sig2.Items));
        [tmp_data2,tmp_ts2,tmp_name2, tmp_unit2] = current_signal(app, tmp);

        tmp_data2(:,2)= tmp_ts2(1):tmp_ts2(1):tmp_ts2(2);

        

        if tmp_ts1(1) < tmp_ts2(1)
            app.settings.lagsig1 = tmp_data1(:,1);
            app.settings.lagts = tmp_data1(:,2);
            app.settings.lagsig2 = interp1(tmp_data2(:,2),tmp_data2(:,1),tmp_data1(:,2));
        else
            app.settings.lagsig2 = tmp_data2(:,1);
            app.settings.lagts = tmp_data2(:,2);
            app.settings.lagsig1 = interp1(tmp_data1(:,2),tmp_data1(:,1),tmp_data2(:,2));
        end        
        
        if ~strcmp(app.popup_lag_int.Value, 'full')       
            tmp = find(vertcat(app.burst_ints.type) ==1);
            brdrs = [];
            for i= 1:length(tmp)
                if strcmp(app.burst_ints(tmp(i)).name, app.popup_lag_int.Value)
                    brdrs = app.burst_ints(tmp(i)).borders;
                end
            end
    
            tmp = find (app.settings.lagts<brdrs(1),1,'last');
            app.settings.lagsig1(1:tmp) = [];
            app.settings.lagsig2(1:tmp) = [];
            app.settings.lagts(1:tmp) = [];
    
            tmp = find (app.settings.lagts>brdrs(2),1,'first');
            app.settings.lagsig1(tmp:end) = [];
            app.settings.lagsig2(tmp:end) = [];
            app.settings.lagts(tmp:end) = [];
        end

        tmp = find(strcmp(app.popup_lag_bgsig.Value, app.popup_lag_bgsig.Items));
        if tmp > 1
            [tmp_data,tmp_ts,~, ~] = current_signal(app, tmp-1);
            tmp_data(:,2)= tmp_ts(1):tmp_ts(1):tmp_ts(2);
            if exist('brdrs')
                tmp = find (tmp_data(:,2)<brdrs(1),1,'last');
                tmp_data(1:tmp,:) = [];
                
                tmp = find (tmp_data(:,2)>brdrs(2),1,'first');
                tmp_data(tmp:end,:) = [];
            end
            app.settings.lagbgsig = tmp_data;
            app.settings.lagbgsig(:,1) = app.settings.lagbgsig(:,1)-min(app.settings.lagbgsig(:,1));
            app.settings.lagbgsig(:,1) = app.settings.lagbgsig(:,1)./max(app.settings.lagbgsig(:,1));
        else
            app.settings.lagbgsig = nan(1,2);
        end
        
        if app.chkbx_lag_ma_sig1.Value
            tmp = floor(round(app.edt_lag_ma_sig1.Value / mean(diff(app.settings.lagts)))/2); %% change window centered on curren idx
            app.settings.lagsig1 = movmean(app.settings.lagsig1, [tmp,tmp]);
        end

        if app.chkbx_lag_ma_sig2.Value
            tmp = floor(round(app.edt_lag_ma_sig2.Value / mean(diff(app.settings.lagts)))/2); %% change window centered on curren idx
            app.settings.lagsig2 = movmean(app.settings.lagsig2, [tmp,tmp]);
        end
        
        app.sldr_lag_start.Limits = [1,length(app.settings.lagts)];
        app.sldr_lag_start.Value = app.sldr_lag_start.Limits(1);
        app.sldr_lag_end.Limits = [1,length(app.settings.lagts)];
        app.sldr_lag_end.Value = app.sldr_lag_end.Limits(2);
end

if plt 
    
    cla(app.ax_lag_sig1)
    cla(app.ax_lag_sig2)
    cla(app.ax_lag_xcorr)
    cla(app.ax_lag_overlay)


    if ~isnan(app.settings.lagbgsig(1,1))
        mn1 = min(app.settings.lagsig1);
        mx1 = max(app.settings.lagsig1);
        mn2 = min(app.settings.lagsig2);
        mx2 = max(app.settings.lagsig2);
        plot (app.ax_lag_sig1,app.settings.lagbgsig(:,2),app.settings.lagbgsig(:,1)*diff([mn1,mx1])+mn1, Color=[0.7 0.7 0.7])
        plot (app.ax_lag_sig2,app.settings.lagbgsig(:,2),app.settings.lagbgsig(:,1)*diff([mn2,mx2])+mn2, Color=[0.7 0.7 0.7])
    else
        plot (app.ax_lag_sig1,nan,nan, Color=[0.7 0.7 0.7])
        plot (app.ax_lag_sig2,nan,nan, Color=[0.7 0.7 0.7])
    end
    hold (app.ax_lag_sig1,'on')
    hold (app.ax_lag_sig2,'on')

    plot (app.ax_lag_sig1,app.settings.lagts,app.settings.lagsig1, Color=[0 0.4470 0.7410])
    plot (app.ax_lag_sig1,app.settings.lagts(1:app.sldr_lag_start.Value),app.settings.lagsig1(1:app.sldr_lag_start.Value), Color=[0.8,0.8,0.8]);
    plot (app.ax_lag_sig1,app.settings.lagts(app.sldr_lag_end.Value:length(app.settings.lagts)),app.settings.lagsig1(app.sldr_lag_end.Value:length(app.settings.lagts)), Color=[0.8,0.8,0.8]);
    hold (app.ax_lag_sig1,'off')
    xlim (app.ax_lag_sig1,[app.settings.lagts(1),app.settings.lagts(end)])

    plot (app.ax_lag_sig2,app.settings.lagts,app.settings.lagsig2, Color=[0.8500 0.3250 0.0980])
    plot (app.ax_lag_sig2,app.settings.lagts(1:app.sldr_lag_start.Value),app.settings.lagsig2(1:app.sldr_lag_start.Value), Color=[0.8,0.8,0.8]);
    plot (app.ax_lag_sig2,app.settings.lagts(app.sldr_lag_end.Value:length(app.settings.lagts)),app.settings.lagsig2(app.sldr_lag_end.Value:length(app.settings.lagts)), Color=[0.8,0.8,0.8]);
    hold (app.ax_lag_sig2,'off')
    xlim (app.ax_lag_sig2,[app.settings.lagts(1),app.settings.lagts(end)])

end

if update_sig

    if ~isnan(app.settings.lagbgsig(1,1))
        mn1 = min(app.settings.lagsig1);
        mx1 = max(app.settings.lagsig1);
        mn2 = min(app.settings.lagsig2);
        mx2 = max(app.settings.lagsig2);
        app.ax_lag_sig1.Children(4).XData = app.settings.lagbgsig(:,2);
        app.ax_lag_sig1.Children(4).YData = app.settings.lagbgsig(:,1)*diff([mn1,mx1])+mn1;
        app.ax_lag_sig2.Children(4).XData = app.settings.lagbgsig(:,2);
        app.ax_lag_sig2.Children(4).YData = app.settings.lagbgsig(:,1)*diff([mn2,mx2])+mn2;
        
    else
        app.ax_lag_sig1.Children(4).XData = nan;
        app.ax_lag_sig1.Children(4).YData = nan;
        app.ax_lag_sig2.Children(4).XData = nan;
        app.ax_lag_sig2.Children(4).YData = nan;
    end

    app.ax_lag_sig1.Children(3).XData =app.settings.lagts;
    app.ax_lag_sig1.Children(3).YData =app.settings.lagsig1;
    xlim (app.ax_lag_sig1, [app.settings.lagts(1), app.settings.lagts(end)])
     
    app.ax_lag_sig2.Children(3).XData =app.settings.lagts;
    app.ax_lag_sig2.Children(3).YData =app.settings.lagsig2;
    xlim (app.ax_lag_sig2, [app.settings.lagts(1), app.settings.lagts(end)])
    drawnow
end

if update_sldr
    tmp = app.sldr_lag_start.Value;
    app.ax_lag_sig1.Children(2).XData =app.settings.lagts(1:tmp);
    app.ax_lag_sig1.Children(2).YData =app.settings.lagsig1(1:tmp);
    app.ax_lag_sig2.Children(2).XData =app.settings.lagts(1:tmp);
    app.ax_lag_sig2.Children(2).YData =app.settings.lagsig2(1:tmp);

    tmp = app.sldr_lag_end.Value;
    app.ax_lag_sig1.Children(1).XData =app.settings.lagts(tmp:end);
    app.ax_lag_sig1.Children(1).YData =app.settings.lagsig1(tmp:end);
    app.ax_lag_sig2.Children(1).XData =app.settings.lagts(tmp:end);
    app.ax_lag_sig2.Children(1).YData =app.settings.lagsig2(tmp:end);
    drawnow
end

if update_ov
    idx = round(app.settings.lag/mean(diff(app.settings.lagts)));
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
    
    [rp, ~,rs, ~] = get_corrs (app,idx+diff(tmp)+1);
    
    if app.settings.is_normal
        if abs(rp - rs) <0.1
            tmp =rp;
        else
            tmp =rs;
        end
    else
        tmp =rs;
    end

    if tmp<0
        sig2 = sig2*(-1);
    end

    plot(app.ax_lag_overlay, sig1, LineWidth=1.5)
    hold (app.ax_lag_overlay, 'on') 
    plot(app.ax_lag_overlay, sig2, LineWidth=1.5)
    hold (app.ax_lag_overlay, 'off')
    xlim (app.ax_lag_overlay, [1,max(length(sig1),length(sig2))])

end
if calc_corr
    tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
    [rp,pp,rs, ps] = get_corrs (app,idx+diff(tmp)+1);
%     tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];

    if app.settings.is_normal
        tmp_str = 'data is normal distributed';
    else
        tmp_str = 'data is NOT normal distributed';
    end
    app.lbl_lag_xcorr.Text = {['Pearson corr: ' num2str(rp,3) ', p: ' num2str(pp,5)],['Spearmans rho: ' num2str(rs,3) ', p: ' num2str(ps,5)],tmp_str};

end

end

function [xc_ts, xc_d, lag] = calc_lac(sig1,sig2,app)
[h1, ~] = kstest((sig1 - mean(sig1)) / std(sig1));
[h2, ~] = kstest((sig2 - mean(sig2)) / std(sig2));

app.settings.is_normal = h1 == 1 && h2 == 1;

fs = 1/mean(diff(app.settings.lagts));
if app.chkbx_lag_wavelet.Value
    waveletName = 'db4'; % Daubechies wavelet
    level = 4; % Decomposition level
    % Perform Wavelet Decomposition
    [C1, L1] = wavedec(sig1, level, waveletName);
    [C2, L2] = wavedec(sig2, level, waveletName);

    % Extract Approximate Coefficients at the coarsest level
    A1 = appcoef(C1, L1, waveletName, level);
    A2 = appcoef(C2, L2, waveletName, level);

    % Compute Cross-Correlation of Approximate Coefficients
    [xcorr_vals, lags] = xcorr(A1, A2, 'coeff');

    % Find the peak in cross-correlation
    [~, max_idx] = max(xcorr_vals);
    
    delay_samples = lags(max_idx);
    % Convert to time delay
    lag = delay_samples / fs;
    xc_ts = lags /fs;
    xc_d = xcorr_vals;
else
    
    [xc_d,xc_ts] = xcorr(sig1, sig2);%abs(xcorr(sig1, sig2));%app.settings.lagxc = abs(xcorr(app.settings.lagsig1(tmp(1):tmp(2)), app.settings.lagsig2(tmp(1):tmp(2))));
%     app.settings.lagxcts = ((-1)*(length(app.settings.lagts(tmp(1):tmp(2)))-1):(length(app.settings.lagts(tmp(1):tmp(2)))-1))*(mean (diff(app.settings.lagts)));
    xc_ts = xc_ts/fs;
    [pks,lcs] = findpeaks(abs(xc_d));
    [pks2,tmp_idx] = sort(pks,1,'descend');
    lcs = lcs(tmp_idx);
    pks =xc_d(lcs);
    rp = nan(size(lcs));
    pp = nan(size(lcs));
    rs = nan(size(lcs));
    ps = nan(size(lcs));
    for i = 1 :length(lcs)
       [rp(i),pp(i),rs(i), ps(i)] = get_corrs (app,lcs(i));
    end
    
    cordif = mean(abs(rp - rs)) <0.1;

    ovrlp = (length(sig1)-abs(lcs-length(sig1)))/length(sig1);
    tmp_idx = find(ovrlp >=0.5);
    lcs = lcs(tmp_idx);
    pks = pks(tmp_idx);
    rp = rp(tmp_idx);
    rs = rs(tmp_idx);
    pp = pp(tmp_idx);
    ps = ps(tmp_idx);
    if app.settings.is_normal
        if cordif
            [~, tmp_idx] = max(abs(rp));
        else
            [~, tmp_idx] = max(abs(rs));
        end
    else
        [~, tmp_idx] = max(abs(rs));
    end
    
    %[~,idx] = max(abs(app.settings.lagxc));
    idx = lcs(tmp_idx);
    if xc_d(idx) < 0
        app.settings.inverse_lagsig = true;
    else
        app.settings.inverse_lagsig = false;
    end
    lag = xc_ts(idx);
end


end

function [rp,pp,rs, ps] = get_corrs (app,idx)
    tmp = [app.sldr_lag_start.Value,app.sldr_lag_end.Value];
    sig1 = app.settings.lagsig1(tmp(1):tmp(2))-mean(app.settings.lagsig1(tmp(1):tmp(2)));
    sig1 = sig1/max(abs(sig1));

    sig2 = app.settings.lagsig2(tmp(1):tmp(2))-mean(app.settings.lagsig2(tmp(1):tmp(2)));
    sig2 = sig2/max(abs(sig2)); 
    idx(2) = idx-length(sig1);
    if idx(2)>= 0

        sidx = [idx(2)+1,length(sig1);...
                1,length(sig2)-idx(2)];
    else

        sidx = [1,length(sig2)-abs(idx(2));...
                abs(idx(2))+1,length(sig2)];
    end
      
    [rp,pp ]  = corr(sig1(sidx(1,1):sidx(1,2)),sig2(sidx(2,1):sidx(2,2)),'Type','Pearson');
    [rs,ps ]  = corr(sig1(sidx(1,1):sidx(1,2)),sig2(sidx(2,1):sidx(2,2)),'Type','Spearman');
    
end

function time_delay = wavelet_delay(sig1, sig2, fs)
    % Estimates the time delay between two signals using Wavelet Transform
    % Inputs:
    %   sig1, sig2  - Input signals (must have the same length)
    %   fs          - Sampling frequency (Hz)
    %   waveletName - Wavelet type (e.g., 'db4', 'sym5')
    %   level       - Number of decomposition levels
    % Output:
    %   time_delay  - Estimated delay in seconds
    waveletName = 'db4'; % Daubechies wavelet
    level = 4; % Decomposition level
    % Perform Wavelet Decomposition
    [C1, L1] = wavedec(sig1, level, waveletName);
    [C2, L2] = wavedec(sig2, level, waveletName);

    % Extract Approximate Coefficients at the coarsest level
    A1 = appcoef(C1, L1, waveletName, level);
    A2 = appcoef(C2, L2, waveletName, level);

    % Compute Cross-Correlation of Approximate Coefficients
    [xcorr_vals, lags] = xcorr(A1, A2, 'coeff');

    % Find the peak in cross-correlation
    [~, max_idx] = max(xcorr_vals);
    figure
plot (lags,xcorr_vals)
    delay_samples = lags(max_idx);

    % Convert to time delay
    time_delay = delay_samples / fs;

    % Display results
    fprintf('Estimated Time Delay: %.6f seconds\n', time_delay);
end