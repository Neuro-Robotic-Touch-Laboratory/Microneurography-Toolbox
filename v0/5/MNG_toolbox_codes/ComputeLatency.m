function [totXrad_mean totXrad_min totXrad_max] = ComputeLatency(radX)

for i=1:length(radX(1,:))
totXrad_mean(i)= mean(nonzeros(radX(:,i)));

if sum(radX(:,i))>0
totXrad_min(i)= min(nonzeros(radX(:,i)));
totXrad_max(i)= max(nonzeros(radX(:,i)));
else
    totXrad_min(i)=0;
    totXrad_max(i)=0;
end

end
totXrad_mean(isnan(totXrad_mean)) = []; 
totXrad_min=nonzeros(totXrad_min)';
totXrad_max=nonzeros(totXrad_max)';

end