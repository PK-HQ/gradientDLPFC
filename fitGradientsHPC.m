function fitGradientsHPC(folderPathEMPlot,...
    folderPathNeuronIdx,referenceElectrodeMap,desiredStatistic,stat,AP,Monkey,uniqueID)
binSize=0.5;
overlapSize=binSize/3;%binSize/2; %dont change
smoothing=1;
% XY limits
minX=-5;
maxX=15;%xticks(-6:2:16);

%end boundaries for each region
a8AD=2.8;
a946D=9.5;
a46D=13.5;
a10D=15.9;
FPd=21.3;

a8AV=3;
a946V=9.7;
a46V=13.8;
a10V=15.9;
FPv=21.3;

a6DR=3.5;
a8B=7.9;
a9=13.9;
a10SUP=14.5;
FPsup=16.4;
%plotIdx={1; 5; 9; 12; 8; 10; 7; 3; 4; 11; 6; 2};
plotIdx={1 2 3 4};
ylimCont=[0 99;...
    60 100;60 100;...
    0 100;0 100;...
    0 100;0 100;0 100;...
    1 8;1 8;...
    0 8;...
    -10 45;-10 45;...
    0 1;...
    50 120;...
    0 1;...
    0 1;...
    0 1;...
    0 1;...
    0 1;...
    0 1;...
    0 1]; %wrong ylim for mixsel

ylimSpacingCont=[33;20;20;...
    20;20;...
    20;20;20;...
    1;1;...
    1;...
    5;5;...
    1;...
    10;...
    1;...
    1;...
    1;...
    1;...
    1;...
    1;...
    1]; %wrong ylim for mixsel

maxVals=[repmat(100,8,1);...
    8;8;...
    100;...
    100;100;...
    1;...
    500];
%% Reference electrode map 632 x 6
%   monkey |  arrayNo |     arrayAB    |  electrodeNo  |     unitNo      |  sessionNo
%    P=1   |   1-12   |  A=1 B=2 nil=0 |      1-32     |  1-5 (useless)  |     1-8
%anatLocations={'vFEF';'dFEF';'v46';'dDLPFC';'8b';'SEF';'vDLPFC';'d46';'a8';'fef';'dDLPFC';'vDLPFC'};
%                1        2     3       4      5    6       7      8    9     10     11        12 
statNames={'All','Tarselpct','Disselpct','D1selpct','D2selpct','CS','LMS','NMSpct',...
    'Trf','D1rf','D2rf','Mixedstr','selF','selFD1','selFD2','selIdx','lat','disGating','lw','selD1','selD2','selD1D2'};
yLab={'All','D1resp','D2resp','Sel','D2sel','CS','LMS','NMS',...
    'Memory field size (D1 & D2)','Memory field size (D1)','Memory field size (D2)','Mixed selectivity',...
    'Selectivity F-stat','Norm. PEV (D1)','Norm. PEV (D2)','Stimulus selectivity index','Response latency (s)',...
    'Distractor filtering','Memory/Motor ratio','D1 selectivity index','D2 selectivity index','D1 & D2 selectivity index'};
xLab={'AP location (mm)'};

fileNamePlots=statNames;

folderNamePlots=statNames;
%                1        2     3       4      5    6       7      8    9     10    11     12       13
%fprintf('Plotting %s\n',statNames{stat})
%merge each electrode with its statistic of interest
%monkeysDataRef=[referenceElectrodeMap,refAllNeurons];
if AP==1
    %monkeysData=[referenceElectrodeMap,desiredStatistic];
    monkeysDataBoth=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,7:8)];
elseif AP==0
    %monkeysData=[referenceElectrodeMap,desiredStatistic];
    monkeysDataBoth=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,8)];
end

%remove duplicate cells
load('data/neuronidx/duplications.mat');
idx=find(dup==1);
monkeysDataBoth(idx,:)=[];

%desiredStatisticBinned=;
%grouping=desiredStatistic;

%ergion-array parcels

dmLPFC=[1 2 4 8 10 11 12]; %[2 4 8]; %[2 4 8 10 11 12];
vmLPFC=[3 7 13]; %[1 3 7]; %[1 3 7 13];
dLPFC=[5 6 9]; %[5 6 9];

region2Fit=1;
%wantedregions={'LPFC','dorso-LPFC','dDLPFC','vDLPFC'};
wantedregions={' ','dDLPFC','vDLPFC'};
wantedregionStrs={'dL','midD','midV'};
ylabel(yLab{stat})
xlabel(xLab{1})
[greenOut,magentaOut]=getColors(9,1);
colors=[greenOut(4,:); greenOut(6,:); greenOut(8,:)];

for reg=[2 3] %dDLPFC and vDLPFC
    wantedregion=wantedregions{reg};fprintf('%s\n',wantedregion);
    monkeysData=monkeysDataBoth(find(monkeysDataBoth(:,1)==Monkey),:);
    switch wantedregion
        case {'dDLPFC'}
            [regionNeurons,~,~]=find(monkeysData(:,end)>=0 & monkeysData(:,end)<8);
            breakpoints=[a8AD,a946D,a46D,a10D];
            allBreaks=[minX,a8AD,a8AD,a946D,a946D,a46D];
            breakpointStrs={'area 8Ad','area 9/46d','area 46d',' area 10md'};
            dataColor=colors(2,:);
            boundarySets=[minX breakpoints(1); breakpoints(1) breakpoints(2); breakpoints(2) breakpoints(3); minX maxX];
        case {'vDLPFC'}
            [regionNeurons,~,~]=find(monkeysData(:,end)<0);
            %[regionNeurons,~,~]=find(monkeysData(:,2)==vmLPFC);
            breakpoints=[a8AV,a946V,a46V,a10V];
            allBreaks=[minX,a8AV,a8AV,a946V,a946V,a46V];
            breakpointStrs={'area 8Av','area 9/46v','area 46v',' area 10mv'};
            dataColor=colors(3,:);
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


    %sliding win smooth
    [smoothedValCont, reductionCount, reductionPercent]=verticalSmoothingGradientSlide(actualData,refData,binSize,overlapSize,stat,colGroup,minX,maxX);
    x=smoothedValCont(2,:)';
    desiredStatistic=smoothedValCont(1,:)';


    %% Generate surrogate data, 3-regions
    %run shuffle data, all 3 regions
    fitDisModelIterHPC_shuffle(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey)
    %run 1seg, all 3 regions
    fitDisModelIterHPC_indiv_rate(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey);
    %run 2seg, all 3 regions
    fitDisModelIterHPC_indiv_rate_2seg(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey);
    %run 3seg, all 3 regions
    fitDisModelIterHPC_indiv_rate_3seg(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey);



    %% Generate surrogate data, 2-region;
    %run shuffle data, 2 regions only
    fitDisModelIterHPC_shuffle_2regions(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey)

    %run 1seg, 2 regions only
    fitDisModelIterHPC_indiv_robust_2regions(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey);
    %run 1seg, 2 regions only    
    fitDisModelIterHPC_indiv_robust_2regions_2seg(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    fileNamePlots{stat},wantedregion,dataColor,wantedregionStrs{reg},folderNamePlots{stat},uniqueID,Monkey);
    
end


