function results = stepwise_regression(app)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ds = app.edt_reg_downsampling.Value;
channel_idx = find(strcmp(app.popup_signal_to_predict.Value,app.popup_signal_to_predict.Items));
[data_p,ts_p,name_p, ~] = current_signal(app, channel_idx);
name_p =char(string(name_p)); 
% data_p_d = downsample(data_p,(.01/ts_p(1)*ds));
[data_p_d, ~] = resample(data_p,ts_p(1):ts_p(1):ts_p(2),100,'linear');

reg_table = app.tbl_regressor.Data;

data_r_d = nan(length(data_p_d),size(reg_table,1));
reg_name = {[]};
lags = nan(length(data_p_d),1);
for i= 1: size(reg_table,1)
    ch_idx = find(strcmp(reg_table{i,1},app.popup_regressor.Items));
    [data_tmp,ts_tmp,name_tmp, unit_tmp] = current_signal(app, ch_idx);
    
%     data_r_d(:,i) = downsample(data_tmp,.01/ts_tmp(1)*ds)';
    [data_r_d(:,i), ~] = resample(data_tmp,ts_tmp(1):ts_tmp(1):ts_tmp(2),100,'linear');
    reg_name{i} = char(string(name_tmp));
    lags(i) = reg_table{i,2};
end

tmp =find(vertcat(app.burst_ints.type) == 1);
borders(1,:) = [1 length(data_p_d)];
int_name{1} = 'full';
for i = 1: length(tmp)        
    borders(i+1,:) = app.burst_ints(tmp(i)).borders;
    borders(i+1,:) = [ceil(borders(i+1,1)/(.01*ds)), floor(borders(i+1,2)/(.01*ds))];
    int_name{i+1} = app.burst_ints(tmp(i)).name;
end

for i = 1: size(borders,1)
    data_p_int = data_p_d(borders(i,1):borders(i,2));
    data_r_int = data_r_d(borders(i,1):borders(i,2),:);

   
    Lag_M=max(lags);
    lag01=Lag_M*1/(.01*ds);

    pred_XXX=double(data_p_int(lag01+1:end));
    pred_XXX(isnan(pred_XXX))=0;
    
    if app.chkbx_reg_mov_average.Value
        pred_XXX=movmean(pred_XXX,app.edt_reg_mov.Value/(.01*ds));
    end
    Regr_XXX= [];
    for j=1:size(data_r_int,2)
        lag02=lags(j)*1/(.01*ds);
        Regr_XXX(:,j)=double(data_r_int(lag01-lag02+1:end-lag02,j));
        Regr_XXX(isnan(Regr_XXX(j,:)),j)=0;
    
        if app.chkbx_reg_mov_average.Value
            Regr_XXX=movmean(Regr_XXX,app.edt_reg_mov.Value/(.01*ds));
        end
    
    end
  
    
    results(i).data = stepwise_GDA_v06(Regr_XXX,reg_name,pred_XXX,app.settings.interval(1,1),name_p);
    results(i).int_name = int_name{i};
    results(i).pred_name = name_p;
    results(i).reg_name = reg_name;

end

end

