function [spec_res, update] = spectral_analysis(app)
%SPECTRAL_ANAL Summary of this function goes here
%   Detailed explanation goes here

channel_idx = find(strcmp(app.popup_signal_spect.Value,app.popup_signal_spect.Items));
[data,ts,~, ~] = current_signal(app, channel_idx);
data(:,2) = ts(1):ts(1):ts(2);

spec_res = struct('x_1_1',{[]}, 'y_1_1',{[]},'lbl_1_1',{[]},...
                   'x_2_1',{[]}, 'y_2_1',{[]},'lbl_2_1',{[]},...
                   'x_3_1',{[]}, 'y_3_1',{[]},'lbl_3_1',{[]},...
                   'x_4_1',{[]}, 'y_4_1',{[]},'lbl_4_1',{[]},...
                   'x_1_2',{[]}, 'y_1_2',{[]},'lbl_1_2',{[]},...
                   'x_2_2',{[]}, 'y_2_2',{[]},'lbl_2_2',{[]},...
                   'x_3_2',{[]}, 'y_3_2',{[]},'lbl_3_2',{[]},...
                   'x_4_2',{[]}, 'y_4_2',{[]},'lbl_4_2',{[]});
ds = [app.edt_ds1.Value, app.edt_ds2.Value];


tmp =find(vertcat(app.burst_ints.type) == 1);
borders = nan(length(tmp)+1,2);

borders(1,:) = [1 size(data,1)];

for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    borders(i+1,:) = [ceil(borders(i+1,1)/ts(1)), floor(borders(i+1,2)/ts(1))];
end

min_dur = min (diff(borders,1,2));

if (min_dur/max(ds))< 160
    ds_max = floor((min_dur)/161);
    warndlg(['Downsampling factor is too high for the selected intervall please choose downsampling factor not higher than : ' num2str(ds_max)]) 
    update = false;
else

    for k =1: size(borders,1)
        data_int = data(borders(k,1):borders(k,2),:);
        for  i = 1:2
        
            a = ds(i);
            fs = (1/ts(1))/a;%1e2*0.05*2000/a;
            data_ds = [downsample(data_int(:,1),a),downsample(data_int(:,2),a)];
            nfft = 1024*16;
               
            
            [spec_res.y_1_1{i,k}, spec_res.x_1_1{i,k}] = pwelch(data_ds(:,1),[],100, nfft, fs);
            spec_res.lbl_1_1{i,k} = ['DS = ' num2str(a)];
           % hold([a6],'on'), plot(a6,f1,pxx,'linewidth',1.5,'DisplayName',['DS=' num2str(a)]);
            
            
            [spec_res.y_2_1{i,k}, spec_res.x_2_1{i,k}] = periodogram(data_ds(:,1),[], nfft, fs);
            spec_res.lbl_2_1{i,k} = ['DS = ' num2str(a)];
            %hold([a7],'on'), plot(a7,f2,pxx1,'linewidth',1.5,'DisplayName',['DS=' num2str(a)]);
    
            %fftfunction
            
            if app.chkbx_windowing.Value
                win = hamming(size(data_ds,1));
                X=fft(data_ds(:,1).*win);
            else
                X=fft(data_ds(:,1));
            end
            y=abs(X);
            spec_res.y_3_1{i,k} = (y.^2 / size(data_ds,1))*a;
            spec_res.x_3_1{i,k} = (0:length(X)-1)*fs/length(X);
            spec_res.lbl_3_1{i,k} = ['DS = ' num2str(a)];
            %hold([a8],'on'), plot(a8,fxx2,pxx2,'linewidth',1.5,'DisplayName',['DS=' num2str(a)]);
           
            
            
            Nw=3;  %effective bandwidth
            [spec_res.y_4_1{i,k}, spec_res.x_4_1{i,k}]=pmtm(data_ds(:,1),Nw, nfft,fs);
            spec_res.lbl_4_1{i,k} = ['DS = ' num2str(a)];
        %     hold([a9],'on'), plot(a9,f4,pxx3,'linewidth',1.5,'DisplayName',['DS=' num2str(a)]);
    
        
            order = [app.edt_o1.Value,app.edt_o2.Value,app.edt_o3.Value];
            for j= 1:3
                save_idx = j+(i-1)*3;
                [spec_res.y_1_2{save_idx,k}, spec_res.x_1_2{save_idx,k}] = pyulear(data_ds(:,1),order(j),nfft,fs);
                spec_res.lbl_1_2{save_idx,k} = ['DS = ' num2str(a) '; O = ' num2str(order(j))];
                %hold([a2],'on'), plot(a2,f5,psd,'linewidth',1.5,'DisplayName',['DS=' num2str(a) '; O=' num2str(order)]);
                
                [spec_res.y_2_2{save_idx,k}, spec_res.x_2_2{save_idx,k}] = pburg(data_ds(:,1), order(j), nfft, fs);
                spec_res.lbl_2_2{save_idx,k} = ['DS = ' num2str(a) '; O = ' num2str(order(j))];
                %hold([a3],'on'), plot(a3,f6,P1,'linewidth',1.5,'DisplayName',['DS=' num2str(a) '; O=' num2str(order)]);
                
                [spec_res.y_3_2{save_idx,k}, spec_res.x_3_2{save_idx,k}] = pcov(data_ds(:,1), order(j), nfft, fs);
                spec_res.lbl_3_2{save_idx,k} = ['DS = ' num2str(a) '; O = ' num2str(order(j))];
                %hold([a4],'on'), plot(a4,f7,P2,'linewidth',1.5,'DisplayName',['DS=' num2str(a) '; O=' num2str(order)]);
                
                [spec_res.y_4_2{save_idx,k}, spec_res.x_4_2{save_idx,k}] = pmcov(data_ds(:,1), order(j), nfft, fs);
                spec_res.lbl_4_2{save_idx,k} = ['DS = ' num2str(a) '; O = ' num2str(order(j))];
                %hold([a5],'on'),
                %plot(a5,f8,P3,'linewidth',1.5,'DisplayName',['DS=' num2str(a) '; O=' num2str(order)];
            
            end 
        end
    end
    update = true;
end
end

