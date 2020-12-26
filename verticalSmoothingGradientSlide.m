function [smoothedValCont, reductionCount, reductionPercent]=verticalSmoothingGradientSlide(actual632,ref632,binSize,overlapSize,stat,colGroup,minX,maxX)
gapSize=binSize-overlapSize; %newly added
minLimX=minX:gapSize:maxX-binSize;
maxLimX=minX+binSize:gapSize:maxX;
valueCol=1;
distanceCol=2;
smoothedValCont=zeros(3,numel(minLimX));
cont=[];
for i=1:numel(minLimX)
    minLim=minLimX(i);
    maxLim=maxLimX(i);
    columnsToSmooth=find(actual632(:,distanceCol)>=minLim & actual632(:,distanceCol)<maxLim);
    regionID=colGroup(find(actual632(:,distanceCol)>=minLim & actual632(:,distanceCol)<maxLim),1);
    cont(i,1:2)=[mean([minLim,maxLim]),numel(columnsToSmooth)];
    if stat<=11
        %for measures that apply to all neurons
        %'D1resp','D2resp','D1sel','D2sel','CS','LMS','NMS','D1rf','D2rf'
        if stat<=8
            %for measures that are proportions
            percentageFactor=100;
        else
            %for measures that are not proportions
            percentageFactor=1;
        end
        if isempty(columnsToSmooth)==1
            numerator=NaN;
        elseif isempty(columnsToSmooth)==0
            numerator=nansum(actual632(columnsToSmooth,valueCol));
        end
        denominator=nansum(ref632(columnsToSmooth,valueCol));
        smoothedVal=numerator*percentageFactor/denominator;
    elseif stat>=12
        %for measures that dont apply to all neurons
        %'Mixedstr','D1pev','D2pev','selIdx','disGating'
        if isempty(columnsToSmooth)==1
            numerator=NaN;
        elseif isempty(columnsToSmooth)==0
            numerator=nanmean(actual632(columnsToSmooth,valueCol));
        end
        smoothedVal=numerator;
        %for measures that are not proportions
        percentageFactor=1;
    end

    
    if ~isempty(unique(regionID))==1
        regionID=unique(regionID);
    else
        regionID=NaN;
    end
        
    %smoothedValCont(1:3,i)=[smoothedVal;mean([minLim maxLim]);regionID];
    smoothedValCont(1:3,i)=[smoothedVal;mean([minLim maxLim]);1];
end

smoothedValCont = smoothedValCont(:,~isnan(smoothedValCont(1,:)));
smoothedValCont=smoothedValCont(:,2:end-1);%remove edges to avoid edge effect
sumValidNow=sum(~isnan((smoothedValCont(1,:))) & (smoothedValCont(1,:) ~= 0));
sumValidPrev=sum(~isnan((actual632(1,:))) & (actual632(1,:) ~= 0));
reductionCount=sumValidNow-sumValidPrev;
reductionPercent=100-((sumValidNow/sumValidPrev)*100);

%{
% plot histogram
figure
histogram(nonzeros(cont(:,2)),'FaceColor',[.5 .5 .5],'LineWidth',2,'BinWidth',1)
%bar(cont(:,1), cont(:,2),'FaceColor',[.5 .5 .5])
%addSkippedTicks([1:1:16],'x')
xticks(0:5:55)
addSkippedTicks([0:1:5],'y') %addSkippedTicks([0:2:12],'y')
xlim([0.48 55.55])
ylim([0 5])
ylabel('Frequency')
xlabel('Neuron count')
upFontSize(22,0.015)
removeURTicks
%}
end