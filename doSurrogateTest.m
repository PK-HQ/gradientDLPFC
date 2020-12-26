function [fitsCell,percentileR]=doSurrogateTest(fitsCell,adjR2Cont,comparisonDistribution,whichSegmentedModel)
% compares actual adjR2 of specified model to the comparison distribution
% e.g. actual 1-seg adjR2 to a dist of 1-segment models. Dist has matched mean and variance 
% of actual data, and with similar adjR2 to the 2-seg model
percentileR=invprctile(comparisonDistribution,adjR2Cont(whichSegmentedModel))/100; %>95%
percentileR=1-percentileR; %inverse to get pvalue <.05
switch percentileR<.05
    case {1}
        fitsCell(whichSegmentedModel)={'+'};
    case {0}
        %no change
        fitsCell(whichSegmentedModel)={'-'};
end