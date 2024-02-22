function hb_res = hb_det_2(varargin)%ecg, threshold, inverse, bp, use_bp)
%HB_DET_2 Summary of this function goes here
%   Detailed explanation goes here

ecg = varargin{1};
threshold = varargin{2};
inverse = varargin{3};
if length(varargin) > 3
    bp = varargin{4};
    use_bp = varargin{5};
    if isempty(bp)
        use_bp = false;
    else
        bp = varargin{4};
    end
end

if use_bp
    ts_ecg = mean(diff(ecg(:,2)));
    ts_bp = bp.ts(1);
    ups= ts_ecg(1)/ts_bp(1);
    ft_idx = bp.foot_idx/ups;
    dcrt_idx = bp.dicrotic_idx/ups;
    if ft_idx(1)< dcrt_idx(1)
        dcrt_idx = [1,dcrt_idx];
    end
    if ft_idx(end)< dcrt_idx(end)
        ft_idx = [ft_idx,size(ecg,1)];
    end
    
end


ecg(:,1) = ecg(:,1)-mean(ecg(:,1));
if inverse
    ecg(:,1) = -ecg(:,1);
end
ecg_dx = -[diff(ecg(:,1));0];
ecg_stdx = std(ecg_dx);

ecg_dx_d = diff(ecg_dx);
tmp_locs = find( ecg_dx_d(:).*circshift(ecg_dx_d(:),[-1,0]) <=0)+1;
locs = tmp_locs(ecg_dx(tmp_locs)>threshold*ecg_stdx);
[pks,locs_2]= findpeaks(ecg_dx, 'MinPeakHeight',threshold*ecg_stdx, 'MinPeakDistance', 200);

if locs_2(1) <121
    locs_2(1) =[];
end
% if locs_2(end) >= size(ecg,1)-20
%     locs_2(end) = [];
% end
ecg_peaks_1st = nan(1,121);

if locs_2(1)-120 < 1
    tmp = ecg(1:locs_2(1));
else
    tmp = ecg(locs_2-120:locs_2(1));
end
ecg_peaks_1st (end-length(tmp)+1:end)= tmp;
for i= -120:0
    ecg_peaks(:,i+121) = ecg(locs_2(2:end)+i,1);
end

ecg_peaks = [ecg_peaks_1st; ecg_peaks];

[~,tmp] = max(ecg_peaks,[],2);
hb_idx = locs_2 +(tmp-120);

if use_bp
    %create indicator colums for test array => odd:dcrt_idx, even: ft_idx 
    test_array = [dcrt_idx',ones(size(dcrt_idx'))*(-1),(1:2:length(dcrt_idx)*2)';
                  hb_idx,zeros(size(hb_idx)),nan(size(hb_idx));
                  ft_idx', ones(size(ft_idx')),(2:2:length(ft_idx)*2)'];
    [~,idx] = sort(test_array(:,1));
    test_array = test_array(idx,:);
    no_bts = find(diff(test_array(:,2)) == 2);
    no_bts_idx_start = test_array(no_bts,1);
    no_bts_idx_end = test_array(no_bts+1,1);
    add_idx = nan(size(no_bts));
    for i = 1:length(no_bts_idx_start)
        [~,lcs,~,p] =  findpeaks(ecg(no_bts_idx_start(i):no_bts_idx_end(i),1));
        [~,mx_idx] = max(p);
        add_idx(i) =  no_bts_idx_start(i)-1+lcs(mx_idx);
    end
    hb_idx = sort([add_idx;hb_idx]);
    mult_bts = find(diff(test_array(:,2)) == 0);
    mult_bts_idx_start = [test_array(mult_bts,1),test_array(mult_bts-1,3)];
    mult_bts(isnan(mult_bts_idx_start(:,2))) = [];
    mult_bts_idx_start(isnan(mult_bts_idx_start(:,2)),:) = [];
    mult_bts(mod(mult_bts_idx_start(:,2),2)==0) = [];
    mult_bts_idx_start(mod(mult_bts_idx_start(:,2),2)==0,:) = [];
    mult_bts_idx_end = test_array(find(ismember(test_array(:,3),mult_bts_idx_start(:,2)+1))-1);
    for i = 1 : length(mult_bts_idx_end)
        tmp_idx = find(hb_idx>=mult_bts_idx_start(i,1) & hb_idx<=mult_bts_idx_end(i));
        [~,mx_idx] = max(ecg(hb_idx(tmp_idx)));
        tmp_idx(mx_idx) = [];
        hb_idx(tmp_idx) = [];
    end
    
    fls_bts = find (diff(test_array(:,2)) == (-1));
    fls_bts(2:2:length(fls_bts))= [];
    fls_bts_idx_start = [test_array(fls_bts+1,1),test_array(fls_bts,3)];
    fls_bts_idx_end = test_array(ismember(test_array(:,3),fls_bts_idx_start(:,2)+1));
    for i = 1 :length(fls_bts_idx_end)
        hb_idx(hb_idx>=fls_bts_idx_start(i,1) & hb_idx<=fls_bts_idx_end(i)) = [];
    end

end

hb_res.idx = hb_idx;
hb_res.t_events = ecg(hb_res.idx,2);


% hb_res.dt_instantaneous = diff(hb_res.t_events);
% hb_res.rate_instantaneous = 1/hb_res.dt_instantaneous * 60;
% hb_res.rate_mean = mean(hb_res.rate_instantaneous);
% hb_res.dtDifferences = diff(hb_res.dt_instantaneous);

%hb_res.rate_meanTot ?????
%hb_res.rateVariability  ????
%hb_res.rateVariability2  ????

%bla = true 




