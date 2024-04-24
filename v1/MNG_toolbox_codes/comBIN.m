function [spikeMatrix, peak_idxs, peakTimes_ms, classes, peak_vals,  dip1, dip2 ]=comBIN(aaa,data,sr,spk_pos,NBINS_peak,NBINS_dip,NBINS_dip2,npointsbeforepeak,npointsafterpeak)
                                                                                  %comBIN(spks,data,sr,1,1,1)
if nargin<8
  npointsbeforepeak=20;
  npointsafterpeak=20;
end

% load(filename,'data','sr');
data=data';

templateMatchedSpikes=aaa.spikes;
peakTimes_ms=aaa.index;
peak_idxs=uint32(aaa.index*sr/1000);

clear index clear spikes



%% Extract window around the spikes

spikeMatrix=zeros(npointsbeforepeak+npointsafterpeak+1,numel(peak_idxs));

for i=1:numel(peak_idxs)
    idx=peak_idxs(i);
    spike=data(idx-npointsbeforepeak:idx+npointsafterpeak);
    spikeMatrix(:,i)=spike;
end
numSpikes=size(spikeMatrix,2);
clear i idx spike

%% find for each spike the 2 dips and delta-t between them (5 features). With respect to zero defined as average of ending and starting point?
peak_vals=data(peak_idxs);
rmv=[];
discardedSpikesMatrix=[];

% tic
% for i=1:numSpikes
%     if spk_pos(i)
%         [minPks,minLocs]=findpeaks(gnegate(spikeMatrix(:,i)));%find all local minima
%         minPks = -minPks;
%     else
%         [minPks,minLocs]=findpeaks(gnegate(-spikeMatrix(:,i)));%find all local minima
%         
%     end
%     
%     if numel(minPks)<2
%         rmv=[rmv i];
%         continue
%     end
%     if spk_pos(i)
%         [~,I]=sort(minPks, 'ascend');
%     else
%         [~,I]=sort(minPks, 'descend');
%     end
%     I=I(1:2); I=sort(I,'ascend'); minLocs2=minLocs(I);%get greatest 2
%     minima2=spikeMatrix(minLocs2,i);
%     dip1(i)=minima2(1);
%     dip2(i)=minima2(2);
%     dip1_idx(i)=minLocs2(1);
%     dip2_idx(i)=minLocs2(2);
%     if minLocs(I(1)) >20
%         fails1 = [fails1,i];
%     end
%     if minLocs(I(2)) <22
%         fails2 = [fails2,i];
%     end
% end
% dip1=dip1';dip2=dip2';dip1_idx=dip1_idx';dip2_idx=dip2_idx';
% toc
%%

pos_idx = find(spk_pos);
neg_idx = find(~spk_pos);
dip1_idx = [];
dip2_idx = [];
dip1 = [];
dip2 = [];
[dip1(pos_idx),dip1_idx(pos_idx)] = min(spikeMatrix(1:20,pos_idx));
[dip1(neg_idx),dip1_idx(neg_idx)] = max(spikeMatrix(1:20,neg_idx));
[dip2(pos_idx),dip2_idx(pos_idx)] = min(spikeMatrix(21:end,pos_idx));
[dip2(neg_idx),dip2_idx(neg_idx)] = max(spikeMatrix(21:end,neg_idx));
dip2_idx = dip2_idx+20;
dip1=dip1';dip2=dip2';dip1_idx=dip1_idx';dip2_idx=dip2_idx';

%%


discardedSpikesMatrix=[discardedSpikesMatrix spikeMatrix(:,rmv)];
spikeMatrix(:,rmv)=[];
peak_idxs(rmv)=[];
peak_vals(rmv)=[];
peakTimes_ms(rmv)=[];
dip1(rmv)=[];
dip1_idx(rmv)=[];
dip2(rmv)=[];
dip2_idx(rmv)=[];


clear I minima2 minPks minLocs minLocs2 i;
% peakRelativeLocation=npointsbeforepeak+1;
peakRelativeLocation=npointsbeforepeak;
delta1=peakRelativeLocation-dip1_idx;
delta2=dip2_idx-peakRelativeLocation;

% discard spikes with some negative delta. 
spikesWithSomeNegativeDelta=unique([find(delta1<0); find(delta2<0)],'sorted');
discardedSpikesMatrix=[discardedSpikesMatrix spikeMatrix(:,spikesWithSomeNegativeDelta)];
spikeMatrix(:,spikesWithSomeNegativeDelta)=[];
peak_vals(spikesWithSomeNegativeDelta)=[];
peak_idxs(spikesWithSomeNegativeDelta)=[];
peakTimes_ms(spikesWithSomeNegativeDelta)=[];
dip1(spikesWithSomeNegativeDelta)=[];
dip1_idx(spikesWithSomeNegativeDelta)=[];
dip2(spikesWithSomeNegativeDelta)=[];
dip2_idx(spikesWithSomeNegativeDelta)=[];
delta1(spikesWithSomeNegativeDelta)=[];
delta2(spikesWithSomeNegativeDelta)=[];

numSpikes=size(spikeMatrix,2);
numDiscarded=size(discardedSpikesMatrix,2);


% Scatterplot the peak, dip1 and dip2 
% figure('WindowState','maximized')
% subplot(2,2,1)
% scatter3(peak_vals, dip1, dip2);
% xlabel("Peak")
% ylabel("Dip 1")
% zlabel("Dip 2")
% subplot(2,2,2)
% scatter(peak_vals, dip1);
% xlabel("Peak")
% ylabel("Dip 1")
% subplot(2,2,3)
% scatter(peak_vals, dip2);
% xlabel("Peak")
% ylabel("Dip 2")
% subplot(2,2,4)
% scatter(dip1, dip2);
% xlabel("Dip 1")
% ylabel("Dip 2")
% sgtitle("Peak and dip values")
% savefig(filename+"_sc1");

% Scatterplot the delta1 and delta2
% figure('WindowState','maximized')
% subplot(2,2,1)
% scatter3(peak_vals, delta1, delta2);
% xlabel("Peak")
% ylabel("Delta 1")
% zlabel("Delta 2")
% subplot(2,2,2)
% scatter(peak_vals, delta1);
% xlabel("Peak")
% ylabel("Delta 1")
% subplot(2,2,3)
% scatter(peak_vals, delta2);
% xlabel("Peak")
% ylabel("Delta 2")
% subplot(2,2,4)
% scatter(delta1, delta2);
% xlabel("Delta 1")
% ylabel("Delta 2")
% sgtitle("Peak and delta values")
% savefig(filename+"_sc2")


%% Group in bins
% NBINS_peak=3;
% NBINS_dip=2;
% NBINS_dip2=1;
NBINS_delta=2;


%Equal frequency
[binMax,edgesMax]=BINfrequency(peak_vals,NBINS_peak);
[binMin1,edgesMin1]=BINfrequency(dip1,NBINS_dip);
[binMin2,edgesMin2]=BINfrequency(dip2,NBINS_dip2);
[binDelta1,edgesDelta1]=BINfrequency(delta1,NBINS_delta);
[binDelta2,edgesDelta2]=BINfrequency(delta2,NBINS_delta);




%% Group spikes according to dip-peak-dip
classes=findgroups(binMax,binMin1,binMin2);

% figure('WindowState','maximized')
% 
% for i=1:max(classes)
%     spikeClass=spikeMatrix(:,classes==i);
%     
%     subplot(5,4,i)
%     hold on
%     plot(spikeClass)
% hold on
% ylim([-20 20])
% 
%     xlabel("Time points (0.1ms)")
%     ylabel("Amplitude")
%      
% end
% title("Spike binning")



end


