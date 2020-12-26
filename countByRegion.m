function [regionCount,sumRegionCount]=countByRegion(regionalNeuronsIdx)
regionCount=[];
for reg=1:size(regionalNeuronsIdx,2)
    regionNeurons=numel(nonzeros(regionalNeuronsIdx(:,reg)));
    regionCount=[regionCount;regionNeurons];
end
sumRegionCount=sum(regionCount);
end