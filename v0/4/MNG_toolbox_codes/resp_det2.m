function results = resp_det2(app,data, thresh)
%RESP_DET Summary of this function goes here
%   Detailed explanation goes here

[b,a] = butter(1,[0.2 1]/50,'bandpass');
flt_data = filtfilt(b,a,data(:,1));
[~, locs] =findpeaks(flt_data,'MinPeakHeight', thresh*std(flt_data), 'MinPeakProminence', .3*std(flt_data));
% figure 
% ax(1)=subplot (3,1,1);
% plot (data(:,1))
% ax(2)=subplot (3,1,2);
% plot(flt_data)
% hold on
% plot(locs,flt_data(locs),'*r') 
% 
% line ([1,50000],[thresh*std(flt_data), thresh*std(flt_data)],'Color', 'r')
% line ([1,50000],[.3*std(flt_data), .3*std(flt_data)],'Color', 'r')
% ax(3) = subplot (3,1,3);
% plot (movmean(flt_data,200))
% linkaxes(ax,'x')
locs = [1;locs;length(flt_data )];
start_idx = nan(length(locs)-1,1);
peak_idx = locs;
for i = 1:length(locs)-1
    [tmp_pk,tmp_lc] = findpeaks(-flt_data(locs(i):locs(i+1),1));
    [~,tmp_idx] = max(tmp_pk);
    if isempty(tmp_idx)
        start_idx(i) = nan;
    else
        start_idx(i) = tmp_lc(tmp_idx) + locs(i)-1;
    end
end

peak_idx([1,length(peak_idx)]) =[];
start_idx(isnan(start_idx)) = [];

test_idx = [start_idx,ones(size(start_idx));peak_idx,zeros(size(peak_idx))];
[~, sort_idx]= sort(test_idx(:,1));
test_idx_s = test_idx(sort_idx,:);
err_idx = find(diff(test_idx_s(:,2))==0);

for i = length(err_idx):-1:1
    if test_idx_s(err_idx(i),2) % two neigboring start idx
        [~,rem_idx] = max(data(test_idx_s([err_idx(i),err_idx(i)+1])));
        test_idx_s(err_idx(i)+rem_idx-1,:) = []; 
    else                    %two neigboring peak dix 
        [~,rem_idx] = min(data(test_idx_s([err_idx(i),err_idx(i)+1])));
        test_idx_s(err_idx(i)+rem_idx-1,:) = [];
    end
end
if ~test_idx_s(1,2)
    test_idx_s(1,:) = [];
end
if ~test_idx_s(end,2)
    test_idx_s(end,:) = [];
end
start_idx = test_idx_s(test_idx_s(:,2)==1,1);
resp_idx =  [start_idx(1:end-1),test_idx_s(test_idx_s(:,2)==0,1),start_idx(2:end)];
results.idx = resp_idx;
results.ts = data(resp_idx(:,2),2);
results.dt_inst = diff(resp_idx(:,2));
results.rate_inst = 60./results.dt_inst;
results.dt_mean = mean(results.dt_inst);
results.rate_mean = mean(results.rate_inst);
end

