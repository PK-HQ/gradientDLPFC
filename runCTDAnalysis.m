function runCTDAnalysis(folderPathNeuronIdx,folderPathCTD,folderPathCTDPlot,...
    datasets, bins_overlap, sessions, trials, m, st, datasets_e, etrials, m_e, st_e,...
    anatLocations, maxbins,testCondition, preLoad, nTarLocations, decodingIters)
%Runs specified number of iterations of cross-temporal decoding for lpfc and individual parcels,
%and plots the individual and compiled decoding heatmaps.

%% Load the functional parcellated neurons per region
fprintf('Loading functionally parcellated neuron index...')
regionalNeuronsIdx=loadMat(folderPathNeuronIdx, 'regionalNeuronsIdxBothMonkeys.mat');

% get count of neurons per region
[countRegions,~]=countByRegion(regionalNeuronsIdx);

% define minimum neurons per region for subsampling of neurons during decoding
minNeurons=min(countRegions);

% testing combinations
testConditions={'Correct' 'target';'Correct' 'distractor';...
    'Error' 'target';'Error' 'distractor'};

%% Define filenames
for cond=testCondition %iterate over all, or choose iteration you want
    
    %define condition for decoding (target x correct/error OR distractor x correct/error)
    testCondC=testConditions{cond,1}; %e.g. correct
    Trial_LabelC=testConditions{cond,2}; %e.g. target
    testCondE=testConditions{cond+2,1}; %e.g. error
    Trial_LabelE=testConditions{cond+2,2}; %e.g. target
    fprintf('Cross-temporal decoding of %s location, %s trials\n',Trial_LabelC,testCondC);
    
    %define folder paths and filename suffixes
    [folderPathCTDCond,fileSuffix]=getFilenameSuffix(testCondC,Trial_LabelC,folderPathCTD);
    [folderPathCTDCondE,fileSuffixE]=getFilenameSuffix(testCondE,Trial_LabelE,folderPathCTD);
    filenamePlots={};

    %% Run CTD for all regions except LPFC, n-iters
    for region=[1 2 4 7 8]%1:numel(anatLocations)
        %randsample to the min number of neurons across regions for decoding
        parcellatedNeuronsIdx=datasample(nonzeros(regionalNeuronsIdx(:,region)),...
            minNeurons,'Replace',false);
        performanceMat=zeros(maxbins,maxbins,decodingIters);
        switch preLoad
            case {1}
                filenamePlot=[folderPathCTDCond sprintf('M1M2_%s_%diters%s.mat',anatLocations{region},1000,fileSuffix)];
                filenamePlots{region}=filenamePlot;

                filenamePlotE=[folderPathCTDCondE sprintf('M1M2_%s_%diters%s.mat',anatLocations{region},1000,fileSuffixE)];
                filenamePlotsE{region}=filenamePlotE;
                1;
            case {0}
                    for iter=1:decodingIters
                        tic
                        fprintf('Decoding Iteration: %d ',iter)
                        if nTarLocations==7
                            performance = runCTD(datasets,m,st,sessions,Trial_LabelC,trials,...
                                bins_overlap,etrials,datasets_e,m_e,st_e,testCondC,parcellatedNeuronsIdx);
                        elseif nTarLocations==4
                            performance = runCTD4(datasets,m,st,sessions,Trial_LabelC,trials,...
                                bins_overlap,etrials,datasets_e,m_e,st_e,testCondC,parcellatedNeuronsIdx);
                        end
                        performanceMat(:,:,iter)=performance;
                        fprintf('(%f mins)\n', (toc/60))
                        %fprintf('\n=========================================\n')
                    end
                    figure()
                    imagesc(flipud(mean(performance,3)),[50 70])
                    %fprintf('\n=========================================\n')
                    filenamePlot=[folderPathCTDCond sprintf('v9_PJparc_%sCTD%dIters%s_R.mat',anatLocations{region},1000,fileSuffix)];
                    save(filenamePlot,'performanceMat')
                    filenamePlots{region}=filenamePlot;
        end
    end


    %% Visualize with heatmap    
    
end
%clf
end
