%% MAIN DIRECTORY
dataPath='/Volumes/Users/PK/Desktop/HPC_PK/';
cd(dataPath)

%% FILEPATHS
%filenames
datasetFilename='datasetsFull632_8R.mat';

%filepaths for inputs/outputs
folderPathRawData='data/raw/';
folderPathNeuronIdx='data/neuronidx/';
folderPathPEV='data/pev/';
folderPathCTD='data/ctd/';
folderPathPCA='data/pca/';
folderPathLW='data/lw/';
folderPathEM='data/em/';
folderPathDistPlot='plots/neuronDist/';
folderPathPEVPlot='plots/pev/';
folderPathCTDPlot='plots/ctd/';
folderPathPCAPlot='plots/pca/';
folderPathEMPlot='plots/em/';
folderPathLWPlot='plots/lw/';

%% SELECT DESIRED ANALYSES
% Calculate functional measures per neuron
executeCalcFuncMeasures=0;

% Get electrode array positions from manual image measurements
executeExtractAP=0;

% Plot 2D map of electrode positons per monkey
execute2DMap=0;

% Gradient analyses
executeGradientFits=1; %fit 1/2/3-segment regression models

% Get functional boundaries, sort neurons by region (for CTD analyses)
executeGetFuncParc=0; 

% Cross-temporal decoding
executeCTD=0; %run crosstemporal decoding, and visualize with heatmaps and boxplots

%% LOAD DATA, SEPARATE CORRECT AND ERROR TRIALS
if exist('datasets','var')==1
    %do nothing if already loaded and extracted in matlab workspace
elseif exist('datasets','var')==0
    saveFiles=0; %save a copy of separated correct and error trials if you need it
    maxbins=60;
    timeBins=1:maxbins;
    [datasets,trials,m,st,bins_overlap,sessions,regional,datasets_e,etrials,m_e,st_e]=...
        fetchCorrectAndErrorTrials(folderPathRawData,datasetFilename,timeBins,saveFiles);
    
    %universal params
    [binsTar,binsDis,binsD1,binsD2]=defEventTimebins(bins_overlap);
    anatLocations={'8Av';'8Ad';'46v';'946d';'8b';'6DR';'946v';'46d'}; %8 region labels, after collapsing DLPFC and FEF labelled neurons
end

%% GET FUNCTIONAL MEASURE PER NEURON
switch executeCalcFuncMeasures
    case {1}
        trialLabel=1; %1=target, %2=distractor, 5=reward
        calcFuncMeasures(folderPathNeuronIdx,folderPathLW,folderPathRawData,...
            datasets,trials,bins_overlap,sessions,datasets_e,etrials,anatLocations,...
            trialLabel,binsD1,binsD2,maxbins);
end

%% Get electrode positions from x- y- and rotation angles obtained from surgery pics
switch executeExtractAP
    case {1}
        [electrodeMapping]=getElectrodePositions(electrodeArrayMeasurements);
        electrodeAPpositions=assignNeuronAPpositions(electrodeMapping, neuronElectrodeIDStrs);
        %getAPDistance_1Sep.m
        %getAPDistance_24Aug
end

%% Plot 2D electrode map
switch execute2DMap
    case {1}
        plot2DMap(folderPathEMPlot,folderPathNeuronIdx)
end

%% FIT 1/2/3-SEGMENT GRADIENT MODELS TO FUNCTIONAL MEASURES, IDENTIFY BEST MODEL
switch executeGradientFits
    case {1}
        % Merge data files from the HPC for preloading
        %mergerShuffled
        %mergerGradient
        AP=1;HPC=0;
        plotGradientFits(folderPathNeuronIdx,folderPathEM,folderPathEMPlot,AP,HPC)
end

%% FUNCTIONAL PARCELLATION OF NEURONS BY REGION
switch executeGetFuncParc
    case {1}
        getFunctionalParc(folderPathNeuronIdx,anatLocations)
end

%% CROSS-TEMPORAL DECODING
switch executeCTD
    case {1}
        % Merge data files from the HPC for preloading
        %mergerCTD
        testCondition=1; %1=target (correct and error), 3=distractor (correct and error)
        preLoad=1; %preload from data run on HPC just for plotting
        nTarLocations=4; %4 (in paper, these have sufficient error trials) or 7
        maxbins=60;
        decodingIters=1;
        runCTDAnalysis(folderPathNeuronIdx,folderPathCTD,folderPathCTDPlot,...
            datasets, bins_overlap,sessions, trials, m, st,datasets_e, etrials, m_e, st_e,...
            anatLocations,maxbins,testCondition,preLoad,nTarLocations,decodingIters)
end



