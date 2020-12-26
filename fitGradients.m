function [dDLPFCCont,vDLPFCCont,dDLPFCfitsPCell,vDLPFCfitsPCell]=fitGradients(...
    folderPathEMPlot,folderPathNeuronIdx,referenceElectrodeMap,desiredStatistic,...
    measure,AP,Monkey,nRegionStr,loadPrefilteredData)
%% FILENAMES, PARAMETERS, INITIALIZING VARS
dataFolder='compiled/'; %compiled = compiled with new data from HPC, or 'archive/' (old stable)
pvalues=[1,1,1];
dDLPFCCont={'','',''};
vDLPFCCont={'','',''};
analysisMode='robust';
maxiters=1000;
preload=1;
%smoothing bin params
binSize=.5;
overlapSize=binSize/3;
smoothing=1;
% XY limits
minX=-5;
maxX=15;
% labels
%{
wantedregions={' ','dDLPFC','vDLPFC'};
wantedregions2={' ','dDLPFC','vDLPFC'};
wantedregionStrsNew={' ','dDLPFC','vDLPFC'};
wantedregionStrsOld={' ','dDLPFC','vDLPFC'};
%}
wantedregions={' ','dDLPFC','vDLPFC'};
wantedregions2={' ','dDLPFC','vDLPFC'};
wantedregionStrsNew={' ','dDLPFC','vDLPFC'};
wantedregionStrsOld={' ','dDLPFC','vDLPFC'};
dDLPFCfitsPCell={'','',''};
vDLPFCfitsPCell={'','',''};
% save file names and regions to analyze
switch nRegionStr
    case {'3'}
        segStrShuffle='';
        segStrSurr='';
        regions=3; %3 regions only in dDLPFC for Monkey 1
    case {'2'}
        segStrShuffle='2region_';
        segStrSurr='2region_surrogate_';
        regions=[2 3]; %2 regions in dDLPFC and vDLPFC for Monkey 1 and 2
end
% save file string
if Monkey==1
    MonkyID='P_';
else
    MonkyID='J_';
end
% analysis mode
switch analysisMode
    case{'linear'}
        robustFitting=0;
    case {'robust'}
        robustFitting=1;
        folderWanted='robust/winLP_share_r';
end

%% ANTERIOR ANATOMICAL BOUNDARIES FOR EACH REGION
a8AD=2.8;
a946D=9.5;
a46D=13.5;
a10D=15.9;
%FPd=21.3;

a8AV=3;
a946V=9.7;
a46V=13.8;
a10V=15.9;
%FPv=21.3;

%[a8AD,a946D,a8AV,a946V]=getFunctionalBounds()

%% XY-AXES LIMITS AND LABELS
monkeysShareAxes=1;
[~,yLimitsdDLFPC,yTickIntervaldDLPFC,yLimitsvDLFPC,yTickIntervalvDLPFC,funcMeasureStrs,yLabels,xLabel]=getYplotparam(Monkey,monkeysShareAxes);

%% FORMAT ELECTRODE/NEURON DATA STRUCTURE FOR ANALYSIS
% Reference electrode map 632 x 6
%   monkey |  arrayNo |     arrayAB    |  electrodeNo  |     unitNo      |  sessionNo
%    P=1   |   1-12   |  A=1 B=2 nil=0 |      1-32     |  1-5 (useless)  |     1-8
%anatLocations={'vFEF';'dFEF';'v46';'dDLPFC';'8b';'SEF';'vDLPFC';'d46';'a8';'fef';'dDLPFC';'vDLPFC'};
%                1        2     3       4      5    6       7      8    9     10     11        12 
if AP==1
    monkeysDataBoth=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,7:8)];
elseif AP==0
    monkeysDataBoth=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,8)];
end

%% remove duplicate cells
load('data/neuronidx/duplications.mat');
idx=find(dup==1);
monkeysDataBoth(idx,:)=[];
fprintf('%s\n',funcMeasureStrs{measure})
for reg=regions %1:numel(wantedregions)%1:numel(wantedregions)
    figure;hold on
    if reg==2
        subplotIdx=1:3;
        ylimCont=yLimitsdDLFPC;
        ylimSpacingCont=yTickIntervaldDLPFC;
        regionStr='dDLPFC';
    elseif reg==3
        subplotIdx=4:6;
        ylimCont=yLimitsvDLFPC;
        ylimSpacingCont=yTickIntervalvDLPFC;
        regionStr='vDLPFC';
    end
    wantedregion=wantedregions{reg};fprintf('%s\n',wantedregion);
    monkeysData=monkeysDataBoth(find(monkeysDataBoth(:,1)==Monkey),:);

    switch wantedregion
        case {' '}
            [regionNeurons,~,~]=find(monkeysData(:,end)>=8); %ylims

                breakpoints=[a6DR,a8B,a9,a10SUP];
                allBreaks=[minX,a6DR,a6DR,a8B,a8B,a9];
                breakpointStrs={'area 6DR',' area 8B','   area 9','area 10d'};
                dataColor=colors(2,:);
                boundarySets=[minX breakpoints(1); breakpoints(1) breakpoints(2); breakpoints(2) breakpoints(3); minX maxX];
        case {'dDLPFC'}
                [redPalette,bluePalette]=getColors(9,1);
                colors=[redPalette(6,:); redPalette(8,:); redPalette(8,:)];
                [regionNeurons,~,~]=find(monkeysData(:,end)>=0 & monkeysData(:,end)<8); %ylims
                breakpoints=[a8AD,a946D,a46D,a10D];
                allBreaks=[minX,a8AD,a8AD,a946D,a946D,a46D];
                breakpointStrs={'area 8Ad','area 9/46d','area 46d',' area 10md'};
                dataColor=colors(2,:);
                boundarySets=[minX breakpoints(1); breakpoints(1) breakpoints(2); breakpoints(2) breakpoints(3); minX maxX];
        case {'vDLPFC'}
                [redPalette,bluePalette]=getColors(9,1);
                colors=[redPalette(6,:); redPalette(8,:); redPalette(8,:)];
                %colors=[magentaOut(6,:); magentaOut(8,:); magentaOut(8,:)];
                [regionNeurons,~,~]=find(monkeysData(:,end)<0); %ylims
                %[regionNeurons,~,~]=find(monkeysData(:,2)==vmLPFC);

                breakpoints=[a8AV,a946V,a46V,a10V];
                allBreaks=[minX,a8AV,a8AV,a946V,a946V,a46V];
                breakpointStrs={'area 8Av','area 9/46v','area 46v',' area 10mv'};
                dataColor=colors(1,:);
                boundarySets=[minX breakpoints(1); breakpoints(1) breakpoints(2); breakpoints(2) breakpoints(3); minX maxX];
    end

    monkeysData1=monkeysData(regionNeurons,:);

    x=[monkeysData1(:,end-1)];
    y=[monkeysData1(:,end)];
    desiredStatistic=monkeysData1(:,7);
    regionsId=[x monkeysData1(:,2)];
    gg=1;
    for i=unique(regionsId(:,2))'
        regionId=find(regionsId(:,2)==i);
        colGroup(regionId,1)=gg;
        gg=gg+1;
    end

    x1=x;
    x2=y;

    actualData=[desiredStatistic x1];
    refData=[ones(size(desiredStatistic)) x1];


    if smoothing==1
        %sliding win smooth
        [smoothedValCont, reductionCount, reductionPercent]=verticalSmoothingGradientSlide(actualData,refData,...
            binSize,overlapSize,measure,colGroup,minX,maxX);
        x=smoothedValCont(2,:)';
        desiredStatistic=smoothedValCont(1,:)';
    elseif smoothing==0
        %collapse electrodes
        meanUniqPosDatas=collapsePosition([x y desiredStatistic],'mean');
        switch measure
            case {4}
                desiredStatistic=meanUniqPosDatas(:,3)*100;
            case {8}
                desiredStatistic=meanUniqPosDatas(:,3)*100;
            otherwise
                desiredStatistic=meanUniqPosDatas(:,3);
        end
        x=meanUniqPosDatas(:,1);
        %y=meanUniqPosDatas(:,2);
    end
        
    % outlier rm
    switch analysisMode
        case {'robust'}
            %no removal, but weighted datapoints during fitting
        case {'median'} %median
            [~,outlierIdx]=rmoutlier(desiredStatistic,folderWanted);
            desiredStatistic(outlierIdx)=[];
            x(outlierIdx)=[];
        case {'mean'} %mean
            [~,outlierIdx]=rmoutlier(desiredStatistic,folderWanted);
            desiredStatistic(outlierIdx)=[];
            x(outlierIdx)=[];
    end
    
    %testing fitting
    switch measure
        case {23}
            xVals=min(x):0.1:max(x);
            yVals=linspace(10,0,numel(xVals));
            yVals=yVals+rand(numel(yVals),1)'*2;
            mdl=fitlm(xVals,yVals);
            desiredStatistic = random(mdl,x);
    end



    %% Load HPC iteration data, run shuffle test
    switch preload
        case {1}
            [shuffleData,surr1Seg,surr2Seg,surr3Seg,...
            surr1SegFilename,surr2SegFilename,surr3SegFilename]...
                =loadShuffledNSurrogateData(dataFolder,segStrShuffle,MonkyID,wantedregionStrsNew,...
                funcMeasureStrs,wantedregionStrsOld,segStrSurr,reg,measure,maxiters,nRegionStr);
            
            % calc shuffled test p-value
            [~,pvalues,~,~]=fitDiscontinuousModel(x,desiredStatistic,allBreaks,yLabels{measure},xLabel{1},ylimCont(measure,:),ylimSpacingCont(measure,:),...
                funcMeasureStrs{measure},wantedregion,dataColor,shuffleData,funcMeasureStrs{measure},wantedregionStrsNew{reg},wantedregions2{reg},subplotIdx,Monkey,robustFitting);
        case {0}
            shuffleData=fitDisModelIter(x,desiredStatistic,allBreaks,yLabels{measure},xLabel{1},ylimCont(measure,:),ylimSpacingCont(measure,:),...
                funcMeasureStrs{measure},wantedregion,dataColor,wantedregionStrsNew{reg},funcMeasureStrs{measure},maxiters);
            shuffleDataFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustShuffle/' segStrShuffle 'r' MonkyID wantedregionStrsNew{reg} '_' funcMeasureStrs{measure} '_'...
                num2str(1000) 'iters.mat'];
            save(shuffleDataFilename,'shuffledAdjR2Cont');
            [~,pvalues,~,~]=fitDiscontinuousModel(x,desiredStatistic,allBreaks,yLabels{measure},xLabel{1},ylimCont(measure,:),ylimSpacingCont(measure,:),...
                funcMeasureStrs{measure},wantedregion,dataColor,shuffleData,funcMeasureStrs{measure},wantedregionStrsNew{reg},wantedregions2{reg},subplotIdx,Monkey);
    end
    
    %% p-value correction
    if sum(~isnan(pvalues))==3
        [~, ~, ~, adjPvals]=fdr_bh(pvalues,.05,'pdep','no');
    elseif sum(~isnan(pvalues))==2
        [~, ~, ~, adjPvals]=fdr_bh(pvalues(1:2),.05,'pdep','no');
        adjPvals=[adjPvals,NaN];
    elseif sum(~isnan(pvalues))==1
        adjPvals=pvalues;
    end
    
    %% LOAD SHUFFLED
    %load(shuffleDataFilename);
    %if strcmp(class(shuffledAdjR2Cont),'cell')
    %    shuffledAdjR2Cont=cell2mat(shuffledAdjR2Cont(:,1:3));
    %end
    
    %% RUN SURROGATE TEST ANA PLOT 2-REGION OR 3-REGION MODELS
    switch nRegionStr
        case {'3'}
            [breakpointCellCont,fitsPCell]=fitDiscontinuousModelAdjP_3model(x,desiredStatistic,...
                allBreaks,yLabels{measure},ylimCont(measure,:),ylimSpacingCont(measure,:),funcMeasureStrs{measure},...
                surr1Seg,surr2Seg,surr3Seg,...
                surr1SegFilename,surr2SegFilename,surr3SegFilename,...
                wantedregionStrsNew{reg},subplotIdx,...
                adjPvals,robustFitting,nRegionStr,loadPrefilteredData);
        case {'2'}
            [breakpointCellCont,fitsPCell]=fitDiscontinuousModelAdjP_2model(x,desiredStatistic,...
                allBreaks,yLabels{measure},ylimCont(measure,:),ylimSpacingCont(measure,:),funcMeasureStrs{measure},...
                surr1Seg,surr2Seg,surr3Seg,...
                surr1SegFilename,surr2SegFilename,surr3SegFilename,...
                wantedregionStrsNew{reg},subplotIdx,...
                adjPvals,robustFitting,nRegionStr,loadPrefilteredData);
    end
    switch reg
        case {3}
            vDLPFCCont=breakpointCellCont;
            vDLPFCfitsPCell=fitsPCell;
        case {2}
            dDLPFCCont=breakpointCellCont;
            dDLPFCfitsPCell=fitsPCell;
    end
    hold off
    upFontSize(22,0.025)
    
    %% SAVE
    %saveFigure(['em/gradient/segmentfits/robustrobust/' dataFolder],...
    %    ['17Dec_fixedR_' regionStr '_' nRegionStr 'regions_' funcMeasureStrs{measure} '_' MonkyID],'')
    CF;
end
end





