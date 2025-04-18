function [threshold, index, par, spikes, spk_pos] = detect_spikes(data, sr, par)
%DETECT_SPIKES Summary of this function goes here
%   Detailed explanation goes here

%default config
init_date = now;

threshold = [];
index = [];
spikes = [];
spk_pos = logical([]);
seglen = par.segments_length*sr*60;
num_seg = ceil (length(data)/(seglen)); %% process 5 min segments
dur_seg = floor(length(data)/num_seg);
idx = nan(num_seg,2);
idx(1,1) = 1; 
for i = 1: num_seg
    idx(i,2) = idx(i,1)+dur_seg-1;
    if i ~= num_seg
        idx(i+1,1) = idx(i,1)+ dur_seg;
    end
end
for n = 1:num_seg  
    x = data(idx(n,1):idx(n,2));
        %<----  Add here extra processing of the signal (x)
    [new_spikes, aux_th, new_index, new_spk_pos]  = amp_detect(x, par);
    index = [index, ((new_index)+idx(n,1)-1)/sr*1000];%/sr*1000 + (idx-1)];
    
    %new_index to ms
    spikes = [spikes; new_spikes];
    threshold = [threshold, aux_th];
    spk_pos = [spk_pos, new_spk_pos];
end

current_par = par;
par = struct;
par = update_parameters(par, current_par, 'detect');
par.detection_date =  datestr(now);
par.channels = 1;
end


function [spikes,thr,index, spk_pos] = amp_detect(x, par)
% Detect spikes with amplitude thresholding. Uses median estimation.
% Detection is done with filters set by fmin_detect and fmax_detect. Spikes
% are stored for sorting using fmin_sort and fmax_sort. This trick can
% eliminate noise in the detection but keeps the spikes shapes for sorting.


sr = par.sr;
w_pre = par.w_pre;
w_post = par.w_post;


if isfield(par,'ref_ms')
    ref = floor(par.ref_ms * par.sr/1000);
else
    ref = par.ref; %for retrocompatibility
end

detect = par.detection;
stdmin = par.stdmin;
stdmax = par.stdmax;

if par.sort_order > 0
    xf = filt_signal(x,par.sort_order,par.sort_fmin,par.sort_fmax,par.sr);
else
    xf = x;
end
if par.detect_order > 0
    xf_detect = filt_signal(x,par.detect_order,par.detect_fmin,par.detect_fmax,par.sr);
else
    xf_detect = x;
end

noise_std_detect = median(abs(xf_detect))/0.6745;
noise_std_sorted = median(abs(xf))/0.6745;
thr = stdmin * noise_std_detect;        %thr for detection is based on detect settings.
thrmax = stdmax * noise_std_sorted;     %thrmax for artifact removal is based on sorted settings.

index = [];
sample_ref = floor(ref/2);
% LOCATE SPIKE TIMES
switch detect
    case 'pos'
        nspk = 0;
        xaux = find(xf_detect(w_pre+2:end-w_post-2-sample_ref) > thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = max((xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
        spk_pos = true(size(index)); 
    case 'neg'
        nspk = 0;
        xaux = find(xf_detect(w_pre+2:end-w_post-2-sample_ref) < -thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = min((xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
        spk_pos = false(size(index)); 
    case 'both'
        nspk = 0;
        xaux = find(abs(xf_detect(w_pre+2:end-w_post-2-sample_ref)) > thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = max(abs(xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
        spk_pos = false(size(index)); 
        idx = find(xf_detect(index)>0);
        spk_pos(idx) = true;
end

% SPIKE STORING (with or without interpolation)
ls = w_pre+w_post;
spikes = zeros(nspk,ls+4);

xf(length(xf)+1:length(xf)+w_post)=0;

for i=1:nspk                          %Eliminates artifacts
    if max(abs( xf(index(i)-w_pre:index(i)+w_post) )) < thrmax
        spikes(i,:)=xf(index(i)-w_pre-1:index(i)+w_post+2);
    end
end
aux = find(spikes(:,w_pre)==0);       %erases indexes that were artifacts
spikes(aux,:)=[];
index(aux)=[];
spk_pos(aux)=[];

switch par.interpolation
    case 'n'
        spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
        spikes(:,1:2)=[];
    case 'y'
        %Does interpolation
        spikes = int_spikes(spikes,par);
end
end

function filtered = filt_signal(x,order,fmin,fmax,sr)
[b,a] = ellip(order,0.1,40,[fmin fmax]*2/sr);
filtered = filtfilt(b, a, x);      
end
