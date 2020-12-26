function [usableFits,lowestValidThreshold]=getUsableFitsPrct(actualmean,actualStd,actualR,surrR,nSegments,ydatapointCol)
%% gets usable fits for the lowest possible threshold difference from the 
%% actual data (same for mean stdev and adjR2) that gives 1000 iters (10-30%)

%param
startThresh=0.01;
stepThresh=.01;
endThresh=0.3;
counter=1;

%init cont
usableFitsCount=cell(length(startThresh:stepThresh:endThresh),3);

for thresholdDiff=startThresh:stepThresh:endThresh
    benchmarkTerm=.2;%abs(thresholdDiff*actualR);
    %% filters for fits that 1) have similar adj R2 (X%)
    %adjR2_Xprct=find(abs(cell2mat(surrR(:,nSegments))-actualR)<abs(thresholdDiff*actualR);

    %% filter for 2) similar mean and var of Y-values (X%)
    %calc surrogate mean and var
    surrDataY_Xprct=cell2mat(surrR(:,ydatapointCol)');
    surrmean_Xprct=mean(surrDataY_Xprct);
    surrStd_Xprct=std(surrDataY_Xprct);
    surrR2Delta_Xprct=cell2mat(surrR(:,nSegments))-actualR;
    surrmeanDelta_Xprct=surrmean_Xprct-actualmean;
    surrStdDelta_Xprct=surrStd_Xprct-actualStd;
    usableFits_Xprct=find(abs(surrR2Delta_Xprct')<benchmarkTerm &...
        abs(surrmeanDelta_Xprct)<abs(thresholdDiff*actualmean) &...
        abs(surrStdDelta_Xprct)<abs(thresholdDiff*actualStd));
    usableFitsCount{counter,1}=thresholdDiff;
    usableFitsCount{counter,2}=numel(usableFits_Xprct);
    usableFitsCount{counter,3}=usableFits_Xprct;
    counter=counter+1;
end

%get the lowest possible threshold with >=1000 iterations
lowestValidThresholdIdx=min(find(cell2mat(usableFitsCount(:,2))>=1000));
if isempty(lowestValidThresholdIdx) & size(surrR,1)==1
    %if zero surrogate datasets
    lowestValidThreshold=0;
    usableFits=1;
elseif isempty(lowestValidThresholdIdx) & size(surrR,1)>1 
    %if non-zero surrogate datasets, but <1000 usable fits
    [maxFits,lowestValidThresholdIdx]=max(cell2mat(usableFitsCount(:,2)));
    lowestValidThreshold=usableFitsCount{lowestValidThresholdIdx,1};
    usableFits=usableFitsCount{lowestValidThresholdIdx,3};
    fprintf('Alert insufficient fits, yield = %.2f %%\n',maxFits*100/length(surrR))
else
    %Have 1000 fits with a good  10-30% threshold, save data in container
    lowestValidThreshold=usableFitsCount{lowestValidThresholdIdx,1};
    usableFits=usableFitsCount{lowestValidThresholdIdx,3};
end

    

%if iters > 1000, truncate it to 1000
if length(usableFits)>1000
    usableFits=usableFits(1:1000);
elseif length(usableFits)<1000 %catch error
    %fprintf('N=%.0f\n',length(usableFits))
    %usableFits=1;
elseif isempty(usableFits) %catch error
    usableFits=0;
end
end