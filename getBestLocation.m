function bestLocationForCellCont=getBestLocation(folderPathNeuronIdx,datasets,trials,...
    sessions,datasets_e,etrials,anatLocations,...
    correctOnly,errorOnly,correctAndErrorTrials)
%% Load data
load([folderPathNeuronIdx 'regionalNeuronsIdx.mat']);
bestLocsFilename=[folderPathNeuronIdx 'bestof7Loc7Correct.mat'];
locations=[2 2; 2 3; 2 4; 3 2; 3 4; 4 2; 4 3; 4 4];
nRegions=size(anatLocations);

%mintrials
trialCounts=[];
if correctAndErrorTrials==1
    load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestLocsForCells4Loc_CE.mat')
    conds=1:2;
    locationsToAnalyze=[2 3 5 6];
    minNeurons=20;
    for session=1:8
        for locationCoord=[2 3 5 6]
            sessionData=[trials(session).val];
            %label CT locs
            TLoc=AssignTrialLabel(sessionData,1);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];

            sessionDataE=[etrials(session).val];
            %label CT locs
            TELoc=AssignTrialLabel(sessionDataE,1);
            idxTELoc=find(TELoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTELoc)];

            %label CT locs
            TLoc=AssignTrialLabel(sessionData,2);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];

            sessionDataE=[etrials(session).val];
            %label CT locs
            TELoc=AssignTrialLabel(sessionDataE,2);
            idxTELoc=find(TELoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTELoc)];

        end
    end
elseif correctAndErrorTrials==0 & correctOnly==1
    load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestLocsForCells7Loc_C.mat')
    conds=2;
    locationsToAnalyze=1:7;
    minNeurons=20;
    for session=1:8
        for locationCoord=1:7
            sessionData=[trials(session).val];
            %label CT locs
            TLoc=AssignTrialLabel(sessionData,1);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];

            %label CT locs
            TLoc=AssignTrialLabel(sessionData,2);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];
        end
    end
elseif correctAndErrorTrials==0 & errorOnly==1
    load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestLocsForCells4Loc_E.mat')
    conds=1;
    locationsToAnalyze=[2 3 5 6];
    minNeurons=20;
    for session=1:8
        for locationCoord=[2 3 5 6]
            sessionData=[etrials(session).val];
            %label CT locs
            TLoc=AssignTrialLabel(sessionData,1);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];

            %label CT locs
            TLoc=AssignTrialLabel(sessionData,2);
            idxTLoc=find(TLoc==locationCoord);
            trialCounts=[trialCounts,numel(idxTLoc)];
        end
    end
end
minTrials=min([trialCounts]);

tic
bestLocationForCellCont=[];
for cond=conds
    neuronID=1;
    correctTrial=cond-1;
    if exist(bestLocsFilename, 'file')==2
        load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestof7Loc7Correct.mat','bestLocationForCellCont')
    else
  
        for i=1:nRegions
            sprintf('%s',anatLocations{i})

            %get neurons for region X
            regionNeurons=nonzeros(regionalNeuronsIdx(:,i));

            for session=1:8
                %for each session, extract neuron spiking data across 89 windows for
                %all locations
                neurons=find(sessions(:,1)==session);
                firstNeuron=min(neurons);lastNeuron=max(neurons);
                sessionNeurons=firstNeuron:lastNeuron;

                %get region's neurons in this session
                regionNeuronsInSession=intersect(regionNeurons,sessionNeurons);
                if isempty(regionNeuronsInSession)==1
                    continue
                else
                    for neuron=regionNeuronsInSession'
                        dataC=datasets(neuron,:,:);
                        dataE=datasets_e(neuron,:,:);

                        %% get best location for cell

                        allLocsCont=NaN(1,7);
                        for locationCoord=locationsToAnalyze
                            %get trials for wanted location 
                            [FRTargetLocs]=fetchFiringRateTarLoc(...
                                neuron,trials,etrials,session,locations,locationCoord,dataC,dataE,correctTrial,correctAndErrorTrials,minTrials);
                            %calc best location  1 cell x 1 loc x 5 bins
                            [maxOfBin]=findBestTarLocForCell(...
                                FRTargetLocs);
                            %store across locs
                            allLocsCont(1,locationCoord,1)=maxOfBin;
                        end
                        %findpeakInZ axis, return col (allLocsCont)
                        [maxValue,valueNumber] = max(allLocsCont(:));
                        if maxValue<=0
                            bestLocationForCell=NaN;
                        elseif maxValue>0
                            [neuron,bestLocationForCell,bin] = ind2sub(size(allLocsCont),valueNumber);
                        end
                        %compile
                        bestLocationForCellCont=[bestLocationForCellCont,bestLocationForCell];

                        neuronID=neuronID+1;
                    end
                end
            end
        end
        %save
        save('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestof7Loc7Correct.mat','bestLocationForCellCont')
    end
end
%load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestLocsForCells7Loc_C.mat') % Consider all 7 locations
%load('/Volumes/Users/PK/Desktop/HPC_PK/data/neuronidx/bestLocsForCells4Loc_CE.mat')
end
    
    
    
    
    