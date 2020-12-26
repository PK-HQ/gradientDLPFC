function [fitsCell,fitsPCell,dist1segRef_usableN,dist1segRef_range,dist2segRef_usableN,dist2segRef_range,...
    dist3segRef_usableN,dist3segRef_range]=whosthebest_v2(...
    actualDataY,pValCont,adjR2Cont,surrogateAdjR2Cont1seg,surrogateAdjR2Cont2seg,surrogateAdjR2Cont3seg,...
    surrogate1segFname,surrogate2segFname,surrogate3segFname,statName,nRegionStr,loadPrefilteredData)
%determine the best model with surrogate stats.
% '-' means n.s., so not compared (model implausible)
% ' ' means model lost (model less plausible)
% '+' means model won (model more plausible)
% if all are ' ' without a '+', then winning model is undetermined, and all are marked '+'.

if isempty(surrogateAdjR2Cont1seg)
    surrogateAdjR2Cont1seg={0,0,0,0,0};
end
if isempty(surrogateAdjR2Cont2seg)
    surrogateAdjR2Cont2seg={0,0,0,0,0};
end
if isempty(surrogateAdjR2Cont3seg)
    surrogateAdjR2Cont3seg={0,0,0,0,0};
end
fitsCell_23wins={'ns','ns','ns'}; %surrogate test outcome for 2-seg and 3-seg VS 1-seg
fitsCell_1wins={'ns','ns','ns'}; %surrogate test outcome for 1-seg VS 2-seg and 3-seg
fitsCell={'ns','ns','ns'}; % final result of surrogate test, collated from both cell arrays above
fitsPCell={NaN,NaN,NaN};
usableFits1seg=[];
usableFits2seg=[];
usableFits3seg=[];
%% set filter for adj R2 of fit, and mean & variance of data points here
 %1 if you want data to be filtered, 0 if not (if you loaded pre-filtered data)
if length(surrogateAdjR2Cont1seg)==1000
    filterData=0;
elseif length(surrogateAdjR2Cont1seg)>1000
    filterData=1;
elseif length(surrogateAdjR2Cont1seg)<1000
    filterData=1;
    sprintf('Warning, insufficient surrogate fits!\n')
end
saveFile=1;
pval_1seg=0;
pval_1seg=0;
%.2 for D2sel James
%.2 Tarsel Pancake
%.25 Tarsel James

dist1segRef_usableN=0;
dist1segRef_range=[-1,-1];
dist2segRef_usableN=0;
dist2segRef_range=[-1,-1];
dist3segRef_usableN=0;
dist3segRef_range=[-1,-1];

is1segSignificant=find(pValCont(1)<0.05);

%% 1 & 2 FITS ONLY
switch sum(~isnan(pValCont))
    case {2}
        %% Shuffled p-value test
        switch numel(is1segSignificant) % if 1-seg model is significant in the shuffled test
            case {0}
                %no change to fitsCell
                fitsCell_23wins={'ns','ns','ns'};
            case {1}
                %% 2-seg vs 1-seg
                %Get usable fits for 2V1
                [dist2seg_1seg,usableFits1seg,lowestValidThreshold]=getUsableFits(surrogateAdjR2Cont1seg,adjR2Cont,1,2,actualDataY,filterData,saveFile,surrogate1segFname,nRegionStr,loadPrefilteredData);
                %Get the comparison distribution for 2V1
                %dist2seg_1seg=cell2mat(surrogateAdjR2Cont1seg(usableFits1seg,2));
                % Surrogate test of 2-seg vs 1-seg
                [fitsCell,pval_2seg]=doSurrogateTest(fitsCell,adjR2Cont,dist2seg_1seg,2);

                %% 1-seg vs 2-seg
                %Get usable fits for 1V2
                [dist1seg_2seg,usableFits2seg,lowestValidThreshold]=getUsableFits(surrogateAdjR2Cont2seg,adjR2Cont,2,1,actualDataY,filterData,saveFile,surrogate2segFname,nRegionStr,loadPrefilteredData);
                %Get the comparison distribution for 1V2
                %dist1seg_2seg=cell2mat(surrogateAdjR2Cont2seg(usableFits2seg,1));
                % Surrogate test of 1-seg vs 2-seg
                [fitsCell,pval_1seg]=doSurrogateTest(fitsCell,adjR2Cont,dist1seg_2seg,1);
                
                fitsPCell{1}=pval_1seg;
                fitsPCell{2}=pval_2seg;
                %% compile outcome
                % 2-model outcome
                %fitsCell{1}=fitsCell_1wins{1};
                %fitsCell{2}=fitsCell_23wins{2};
                
                
                %% occasionally used for plotting
                % Get number of usables, range 1seg
                dist1segRef_usableN=numel(usableFits1seg);
                dist1segRef_range(1)=min(cell2mat(surrogateAdjR2Cont1seg(:,1)));
                dist1segRef_range(2)=max(cell2mat(surrogateAdjR2Cont1seg(:,1)));
                
                %get number of usables, range 2seg
                dist2segRef_usableN=numel(usableFits2seg);
                dist2segRef_range(1)=min(cell2mat(surrogateAdjR2Cont2seg(:,2)));
                dist2segRef_range(2)=max(cell2mat(surrogateAdjR2Cont2seg(:,2)));
                
        
        end



    %% 1 & 2 & 3 FITS
    case {3}
        %% Surrogate test
        switch numel(is1segSignificant) % if 1-seg model is significant in the shuffled test
            case {0}
                %no change to fitsCell
                fitsCell_23wins={'ns','ns','ns'};
            case {1}
                %% 2-seg and 3-seg VS 1-seg
                %Get usable fits for 2V1 and 3v1
                [dist2seg_1seg,usableFits1seg,lowestValidThreshold]=getUsableFits(surrogateAdjR2Cont1seg,adjR2Cont,1,2,actualDataY,filterData,saveFile,surrogate1segFname,nRegionStr,loadPrefilteredData);
                [dist3seg_1seg,~,~]=getUsableFits(surrogateAdjR2Cont1seg,adjR2Cont,1,3,actualDataY,filterData,saveFile,surrogate1segFname,nRegionStr,loadPrefilteredData);
                % Get comparison distributions
                %dist2seg_1seg=cell2mat(surrogateAdjR2Cont1seg(usableFits1seg,2));
                %dist3seg_1seg=cell2mat(surrogateAdjR2Cont1seg(usableFits1seg,3));
                % Surrogate test of 2-seg vs 1-seg
                [fitsCell,pval_2v1seg]=doSurrogateTest(fitsCell,adjR2Cont,dist2seg_1seg,2);
                % Surrogate test of 3-seg vs 1-seg
                [fitsCell,pval_3v1seg]=doSurrogateTest(fitsCell,adjR2Cont,dist3seg_1seg,3);
                
                %% 1-seg VS 2-seg and 3-seg 
                %Get usable fits for 1V2
                [dist1seg_2seg,usableFits2seg,lowestValidThreshold]=getUsableFits(surrogateAdjR2Cont2seg,adjR2Cont,2,1,actualDataY,filterData,saveFile,surrogate2segFname,nRegionStr,loadPrefilteredData);
                % Get comparison distribution
                %dist1seg_2seg=cell2mat(surrogateAdjR2Cont2seg(usableFits2seg,1));
                % Surrogate test of 1-seg vs 2-seg
                [fitsCell_1wins,pval_1v2seg]=doSurrogateTest(fitsCell_1wins,adjR2Cont,dist1seg_2seg,2);
                
                %Get usable fits for 1V3
                [dist1seg_3seg,usableFits3seg,lowestValidThreshold]=getUsableFits(surrogateAdjR2Cont3seg,adjR2Cont,3,1,actualDataY,filterData,saveFile,surrogate3segFname,nRegionStr,loadPrefilteredData);
                % Get comparison distribution
                %dist1seg_3seg=cell2mat(surrogateAdjR2Cont3seg(usableFits3seg,1));
                % Surrogate test of 1-seg vs 3-seg
                [fitsCell_1wins,pval_1v3seg]=doSurrogateTest(fitsCell_1wins,adjR2Cont,dist1seg_3seg,3);
                
                switch strcmp(fitsCell_1wins{2},'+') & strcmp(fitsCell_1wins{3},'+')
                    case {1} %if 1 wins 2 and 3, +
                        fitsCell(1)={'+'};
                    case {0} %Else -
                        fitsCell(1)={'-'};
                end
                fitsPCell{1}=[pval_1v2seg,pval_1v3seg];
                fitsPCell{2}=pval_2v1seg;
                fitsPCell{3}=pval_3v1seg;
                %% compile outcome
                % 3-model outcome
                %fitsCell{1}=fitsCell_1wins{1};
                %fitsCell{2}=fitsCell_23wins{2};
                %fitsCell{3}=fitsCell_23wins{3};
                %% occasionally used for plotting
                %get number of usables, range
                dist1segRef_usableN=numel(usableFits1seg);
                dist1segRef_range(1)=min(cell2mat(surrogateAdjR2Cont1seg(:,1)));
                dist1segRef_range(2)=max(cell2mat(surrogateAdjR2Cont1seg(:,1)));            

                %get number of usables, range
                dist2segRef_usableN=numel(usableFits2seg);
                dist2segRef_range(1)=min(cell2mat(surrogateAdjR2Cont2seg(:,2)));
                dist2segRef_range(2)=max(cell2mat(surrogateAdjR2Cont2seg(:,2)));
                
                %get number of usables, range
                dist3segRef_usableN=numel(usableFits3seg);
                dist3segRef_range(1)=min(cell2mat(surrogateAdjR2Cont3seg(:,3)));
                dist3segRef_range(2)=max(cell2mat(surrogateAdjR2Cont3seg(:,3)));
        end
end

end





