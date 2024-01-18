function [totXrad_mean totXrad_min totXrad_max] = ComputeLatencyCell(radX)

for i=1:length(radX)
totXrad_mean(i)= mean(radX{i});

if sum(radX{i})>0
totXrad_min(i)= radX{i}(1);
totXrad_max(i)= radX{i}(end);
else
    totXrad_min(i)=0;
    totXrad_max(i)=0;
end

end
totXrad_mean(isnan(totXrad_mean)) = []; 
totXrad_min=nonzeros(totXrad_min)';
totXrad_max=nonzeros(totXrad_max)';

end