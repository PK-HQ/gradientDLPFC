function [break1,break2,allBreak1,allBreak2]=extract2Breaks(DLPFCstatsCont,DLPFCstatsUsed,statNames,regionStr,MonkyID)
% StatsCont is a cell, with measure x n-segment fit (1, 2 and 3 segments). Function 
% here extracts breakpoints from significant 2seg and 3seg fits, for plotting
% and getting functional boundaries.

oneBreakCont=[];
twoBreakCont=[];

for stat=DLPFCstatsUsed%setdiff(9:22,[9 13 22])%[10 11 12 15 16 17 18]
    twoSegFits=DLPFCstatsCont{stat,2}; %2 break model
    threeSegFits=DLPFCstatsCont{stat,3}; %3 break model
    oneBreakCont(stat,1:2)=[0,0];
    twoBreakCont(stat,1:3)=[0,0,0];
    %% if empty, set all to NaN
    switch isempty(threeSegFits)
        case {1}
            threeSegFits={NaN,NaN,NaN};
    end
    switch isempty(twoSegFits)
        case {1}
            twoSegFits={NaN,NaN,NaN};
    end
    
    %% check if fit is significant. If yes, take the breakpoints
    switch twoSegFits{1}==1
        case {1}
            %get break position
            oneBreakpos=twoSegFits{2};
            oneBreakCont(stat,1:2)=[oneBreakpos,stat]; % store this
        case {0}
            oneBreakpos=0;
    end
    switch threeSegFits{1}==1
        case {1}
            twoBreakpos=threeSegFits{2};
            twoBreakCont(stat,1:3)=[twoBreakpos,stat]; % store this
        case {0}
            twoBreakpos=[0,0];
    end
end

%Get measures with >0 points per line segment
allBreak1=[oneBreakCont(find(oneBreakCont(:,1)<5 & oneBreakCont(:,1)>0),:);...
    twoBreakCont(find(twoBreakCont(:,1)>0),[1 3])]; %get 2seg gradient fits, 1st break
allBreak2=[oneBreakCont(find(oneBreakCont(:,1)>5 & oneBreakCont(:,1)>0),:);...
    twoBreakCont(find(twoBreakCont(:,2)>0),[2 3])]; %get 2seg gradient fits, 2nd break

% find 2 values for same measure
[~, ind] = unique(allBreak1(:,2), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(allBreak1, 1), ind);
if numel(duplicate_ind)>0
    % duplicate values
    duplicateVal = allBreak1(duplicate_ind, 2);
    for i=1:numel(duplicateVal)
        dupe=duplicateVal(i);
        x=allBreak1(find(allBreak1(:,2)==dupe),1);
        meanVal=nanmean(x);
        %replace
        allBreak1(find(allBreak1(:,2)==dupe),:)=[];
        allBreak1=[allBreak1;meanVal,dupe];
    end
end
% find 2 values for same measure
[~, ind] = unique(allBreak2(:,2), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(allBreak2, 1), ind);
if numel(duplicate_ind)>0
    % duplicate values
    duplicateVal = allBreak2(duplicate_ind, 2);
    for i=1:numel(duplicateVal)
        dupe=duplicateVal(i);
        x=allBreak2(find(allBreak2(:,2)==duplicateVal),1);
        meanVal=nanmean(x);
        %replace
        allBreak2(find(allBreak2(:,2)==dupe),:)=[];
        allBreak2=[allBreak2;meanVal,duplicateVal];
    end
end

%Extract mean 2 breaks
break1=nanmedian(allBreak1(:,1));
break2=nanmedian(allBreak2(:,1));
fprintf([MonkyID ' ' regionStr ' 2 Breaks: %.2f (%.2f) and %.2f (%.2f)\n'],break1, break1-2.8,break2,break2-9.5)

end