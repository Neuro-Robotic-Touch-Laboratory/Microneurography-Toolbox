function [radX_c SpikeTimeX_c,first,last] = computeRadFiring(t_ECG,t_SPIKES)
% t_ECG=t_events_plus,
% t_SPIKES,
% ind_SPIKES=1:length(t_SPIKES);


t_CORTOpeak=t_SPIKES';%ECG_RESULTS.t_events;  %BELT_Rpeak;

t_LUNGOpeak=t_ECG;%t_maxima(1:end);

%RR interval 
rrint= diff(t_CORTOpeak);
tint=t_CORTOpeak(2:end);

tmat_c = {[]};
rrmat_c ={[]};
first = nan(length(t_ECG)-1,1);
last = nan(length(t_ECG)-1,1);
for i =1:length(t_LUNGOpeak)-1
    tmat_c{i} = tint(tint>=t_LUNGOpeak(i) & tint<t_LUNGOpeak(i+1));
    rrmat_c{i} = rrint(tint>=t_LUNGOpeak(i) & tint<t_LUNGOpeak(i+1));
    if~isempty(tmat_c{i})
        first(i,1) = tmat_c{i}(1);
        last(i,1) = tmat_c{i}(end);
    else
        first(i,1) = nan;
        last(i,1) = nan;
    end
end

radX_c = {[]};
SpikeTimeX_c = {[]};
DTresp_v = diff(t_LUNGOpeak);
for i = 1 : length(tmat_c)
    radX_c{i} = (tmat_c{i}-t_LUNGOpeak(i))/DTresp_v(i)*2*pi;
    SpikeTimeX_c{i} = tmat_c{i}-t_LUNGOpeak(i);
end

end