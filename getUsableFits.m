function [comparisonDist,usableFits,lowestValidThreshold]=getUsableFits(surrogateDataR,actualR,nSegments,wantedSegment,actualDataY,filterData,saveFile,saveFileName,nRegionStr,loadPrefilteredData)
%% filters for fits that 1) have similar adj R2 AND 2) similar mean and var of Y-values (within 10%-30%)
%% dont filter if pre-filtered!
ydatapointCol=5;
usableFits=[];
lowestValidThreshold=NaN;
switch loadPrefilteredData
    case {1}
        comparisonDist=cell2mat(surrogateDataR(:,wantedSegment));
    case {0}
        switch filterData
            case {1}
                if isempty(surrogateDataR)
                    usableFits=[];
                    lowestValidThreshold=NaN;
                else
                    %% find mean,var,adjR2 that gives 1000 iterations, ideally 10% for all (else, allow % threshold to vary consistently across all 3)

                    %calc actual mean and var
                    actualmean=mean(actualDataY);
                    actualStd=std(actualDataY);
                    actualR=actualR(nSegments);
                    [usableFits,lowestValidThreshold]=getUsableFitsPrct(actualmean,actualStd,actualR,surrogateDataR,nSegments,ydatapointCol);
                    switch saveFile
                        case {1}
                            surrogateAdjR2Cont=surrogateDataR(usableFits,:);
                            save([saveFileName(1:end-4) '_filtered' num2str(lowestValidThreshold*100) 'prct_' num2str(size(surrogateAdjR2Cont,1)) 'iters_' nRegionStr 'seg_withThreshold.2.mat'],'surrogateAdjR2Cont');
                            save([saveFileName(1:end-4) '_filtered_' nRegionStr 'seg.2.mat'],'surrogateAdjR2Cont');
                    end
                end
            case {0}
                %dont filter usable fits
                maxCount=min(size(surrogateDataR,1),1000);
                usableFits=1:maxCount;
                lowestValidThreshold=NaN;
        end

        comparisonDist=cell2mat(surrogateDataR(usableFits,wantedSegment));
end
end