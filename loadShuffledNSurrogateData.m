function [shuffledAdjR2Cont,surrogate1segData,surrogate2segData,surrogate3segData,...
    surrogate1segLoadFilename,surrogate2segLoadFilename,surrogate3segLoadFilename]=loadShuffledNSurrogateData(dataFolder,segStrShuffle,MonkyID,wantedregionStrsNew,folderNamePlots,wantedregionStrsOld,...
    segStrSurr,reg,stat,maxiters,nRegionStr)
wantedregionStrsNew2={'xLPFC','dDLPFC','vDLPFC'};

%% Shuffled distribution filenames
%{
if reg==2
    %filename changes to 2 region shuffled due to segStr2
    shuffleDataFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustShuffle/' segStrShuffle 'r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
        num2str(1000) 'iters.mat'];
    shuffledSaveFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustShuffle/' segStrShuffle 'r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
        num2str(1000) 'iters.mat'];
elseif reg==3
    %no change for vDLPFC
    shuffleDataFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustShuffle/r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
        num2str(1000) 'iters.mat'];
    shuffledSaveFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustShuffle/r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
        num2str(1000) 'iters.mat'];
end
%}
shuffleDataFilenameOld=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/' wantedregionStrsOld{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters.mat'];

%% Surrogate distribution filenames
%{
%Unfiltered and untruncated, to process before analysis
surrogate1segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters.mat'];
surrogate2segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '2seg_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters.mat'];
surrogate3segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '3seg_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters.mat'];
%Filtered and truncated, ready for analysis
surrogate1segLoadFilename_truncated=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered_' nRegionStr 'seg.2.mat'];
surrogate2segLoadFilename_truncated=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '2seg_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered_' nRegionStr 'seg.2.mat'];
surrogate3segLoadFilename_truncated=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/robust/robustSurrogateNew/'...
    dataFolder segStrSurr '3seg_r' MonkyID wantedregionStrsNew{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered_' nRegionStr 'seg.2.mat'];
%}

%Filtered and truncated, ready for analysis
surrogate1segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/surrogateData/' ...
    '2region_' 'r' MonkyID wantedregionStrsNew2{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered.mat'];
surrogate2segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/surrogateData/'...
    '2region_' '2seg_r' MonkyID wantedregionStrsNew2{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered.mat'];
surrogate3segLoadFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/surrogateData/'...
    '2region_' '3seg_r' MonkyID wantedregionStrsNew2{reg} '_' folderNamePlots{stat} '_'...
    num2str(maxiters) 'iters_filtered.mat'];
shuffleDataFilename=['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffledData/' segStrShuffle 'r' MonkyID wantedregionStrsNew2{reg} '_' folderNamePlots{stat} '_'...
        num2str(1000) 'iters_filtered.mat'];
surrogateAdjR2Cont={0,0,0,0,0};

if exist(shuffleDataFilename, 'file') == 2
    % File exists.
    load(shuffleDataFilename);
    if strcmp(class(shuffledAdjR2Cont),'cell')
        shuffledAdjR2Cont=cell2mat(shuffledAdjR2Cont(:,1:3));
        shuffledAdjR2Cont(:,3)=zeros(1000,1);
    end
    shuffledAdjR2Cont=shuffledAdjR2Cont(1:1000,:);
    %save(shuffleDataFilename_truncated_new,'shuffledAdjR2Cont')
elseif exist(shuffleDataFilenameOld, 'file') == 2 % if new file name doesn't exist
    load(shuffleDataFilenameOld);
    fprintf('Shuffled data file doesn''t exist, old only\n')
    %[~,pvalues,~,~]=fitDiscontinuousModel(x,desiredStatistic,allBreaks,yLab{stat},xLab{1},ylimCont(stat,:),ylimSpacingCont(stat,:),...
    %    fileNamePlots{stat},wantedregion,dataColor,shuffledAdjR2Cont,folderNamePlots{stat},wantedregionStrsNew{reg},wantedregions2{reg},subplotIdx,subpltLabel,Monky,robust);
    %save(shuffleDataFilename,'shuffledAdjR2Cont');
else
    fprintf('Shuffled data file doesn''t exist, run fitDisModelIter()\n')
end

% Load pre-filtered 1 seg surrogate data, else load raw
if exist(surrogate1segLoadFilename, 'file') == 2
    load(surrogate1segLoadFilename);
    surrogate1segData=surrogateAdjR2Cont;
else
    surrogate1segData={0,0,0,0,0};
end

% Load pre-filtered 2 seg surrogate data, else load raw
if exist(surrogate2segLoadFilename, 'file') == 2
    load(surrogate2segLoadFilename);
    surrogate2segData=surrogateAdjR2Cont;
else
    surrogate2segData={0,0,0,0,0};
end

% Load pre-filtered 3 seg surrogate data, else load raw
if exist(surrogate3segLoadFilename,'file')==2
    load(surrogate3segLoadFilename);
    if isempty(surrogateAdjR2Cont)
        surrogate3segData={0,0,0,0,0};
    else
        surrogate3segData=surrogateAdjR2Cont;
    end
else
    surrogate3segData={0,0,0,0,0};
end
end