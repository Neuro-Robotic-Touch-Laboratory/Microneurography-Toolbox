function [binIdx,binEdges]=BINfrequency(input_vect,nBins)
%Limitation: it may in separate in different bins repeated numbers. Shouldnt be much problem with real numbers were the p
%of 2 numbers being equal is small.

Q=fix(numel(input_vect)/nBins);
R=rem(numel(input_vect),nBins);
binCounts=[Q*ones(1,R)+1 Q*ones(1,nBins-R)];%nbins with equal number of elements and one extra elements to the 1st R bins if the division is not integer

[~,sortedIdx]=sort(input_vect,'ascend');

binEdges=[];
binIdx=zeros(numel(input_vect),1,'int8')';
first=0;
last=0;
for i=1:numel(binCounts)
    first=last+1;
    last=first+binCounts(i)-1;
    binIdx(sortedIdx(first:last))=i;
    if first==1
        binEdges=[binEdges (input_vect(sortedIdx(first)))];
    else
        binEdges=[binEdges (input_vect(sortedIdx(first-1))+input_vect(sortedIdx(first)))/2];
    end
end
binEdges=[binEdges input_vect(sortedIdx(end))];
end