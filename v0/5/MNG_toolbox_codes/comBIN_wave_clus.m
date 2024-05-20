function [spikeMatrix, peak_idxs, peakTimes_ms, classes, peak_vals,  dip1, dip2 ]=comBIN_wave_clus(cluster,data,sr)

npointsbeforepeak=20;
npointsafterpeak=20;

peakTimes_ms=cluster.cluster_class(:,2)';
peak_idxs=uint32(cluster.cluster_class(:,2)'*sr/1000);

clear index clear spikes
classes = (cluster.cluster_class(:,1))';
%% Extract window around the spikes

spikeMatrix=zeros(npointsbeforepeak+npointsafterpeak+1,numel(peak_idxs));

for i=1:numel(peak_idxs)
    idx=peak_idxs(i);
    spike=data(idx-npointsbeforepeak:idx+npointsafterpeak);
    spikeMatrix(:,i)=spike;
end
numSpikes=size(spikeMatrix,2);
clear i idx spike
%% necessary ??
%% find for each spike the 2 dips and delta-t between them (5 features). With respect to zero defined as average of ending and starting point?
peak_vals=data(peak_idxs);
rmv=[];
discardedSpikesMatrix=[];
for i=1:numSpikes
    [minPks,minLocs]=findpeaks(gnegate(-spikeMatrix(:,i)));%find all local minima
    if numel(minPks)<2
        rmv=[rmv i];
        continue
    end
    [~,I]=sort(minPks, 'descend');I=I(1:2); I=sort(I,'ascend'); minLocs2=minLocs(I);%get greatest 2
    minima2=spikeMatrix(minLocs2,i);
    dip1(i)=minima2(1);
    dip2(i)=minima2(2);
    dip1_idx(i)=minLocs2(1);
    dip2_idx(i)=minLocs2(2);
end
dip1=dip1';dip2=dip2';dip1_idx=dip1_idx';dip2_idx=dip2_idx';

% discardedSpikesMatrix=[discardedSpikesMatrix spikeMatrix(:,rmv)];
% spikeMatrix(:,rmv)=[];
% peak_idxs(rmv)=[];
% peak_vals(rmv)=[];
% peakTimes_ms(rmv)=[];
% dip1(rmv)=[];
% dip1_idx(rmv)=[];
% dip2(rmv)=[];
% dip2_idx(rmv)=[];
% classes(rmv)=[];
% 
% clear I minima2 minPks minLocs minLocs2 i;
% % peakRelativeLocation=npointsbeforepeak+1;
% peakRelativeLocation=npointsbeforepeak;
% delta1=peakRelativeLocation-dip1_idx;
% delta2=dip2_idx-peakRelativeLocation;
% 
% % discard spikes with some negative delta. 
% spikesWithSomeNegativeDelta=unique([find(delta1<0); find(delta2<0)],'sorted');
% discardedSpikesMatrix=[discardedSpikesMatrix spikeMatrix(:,spikesWithSomeNegativeDelta)];
% spikeMatrix(:,spikesWithSomeNegativeDelta)=[];
% peak_vals(spikesWithSomeNegativeDelta)=[];
% peak_idxs(spikesWithSomeNegativeDelta)=[];
% peakTimes_ms(spikesWithSomeNegativeDelta)=[];
% dip1(spikesWithSomeNegativeDelta)=[];
% dip1_idx(spikesWithSomeNegativeDelta)=[];
% dip2(spikesWithSomeNegativeDelta)=[];
% dip2_idx(spikesWithSomeNegativeDelta)=[];
% delta1(spikesWithSomeNegativeDelta)=[];
% delta2(spikesWithSomeNegativeDelta)=[];
% classes(spikesWithSomeNegativeDelta)=[];
% 
% numSpikes=size(spikeMatrix,2);
% numDiscarded=size(discardedSpikesMatrix,2);

%% Group in bins
% classes = (wc_res.cluster_class(:,1))';
% NBINS_delta=2;
% 
% %Equal frequency
% [binMax,edgesMax]=BINfrequency(peak_vals,NBINS_peak);
% [binMin1,edgesMin1]=BINfrequency(dip1,NBINS_dip);
% [binMin2,edgesMin2]=BINfrequency(dip2,NBINS_dip2);
% [binDelta1,edgesDelta1]=BINfrequency(delta1,NBINS_delta);
% [binDelta2,edgesDelta2]=BINfrequency(delta2,NBINS_delta);
% 
% %% Group spikes according to dip-peak-dip
% classes=findgroups(binMax,binMin1,binMin2);
end


