function plotGradientFits(folderPathNeuronIdx,folderPathEM,folderPathEMPlot,AP,HPC)
%% FITS 2-REGION AND 3-REGION FITS TO SMOOTHED DATA ALONG THE AP AXIS, 
%% AND TESTS DISTRIBUTION FITS OF THE MODEL COMPARISON P-VALUES  
wantedRegions={'2','3'};
wantedMonkeys=[1 2];
loadFilteredData=0;

%% Load single neuron data
load([folderPathNeuronIdx 'functionalMeasures.mat']);
funcMeasures=funcMeasures(:,1);
load([folderPathNeuronIdx 'electrodeMappingCurved.mat']);

%% Data structure
%{
1 1:632 neurons 
2 target selective
3 distractor selective
4 delay 1 selective 
5 delay 2 selective 
6x classical selective 
7x linear mixed selective 
8 nonlinear mixed selective

9 Receptive field size 
10 D1 field size 
11 D2 field size 
12x NMS interaction term 
13x NMS main effect 

14x Delay 1 selective f-measure
15x Delay 2 selective f-measure 
16 Selectivity index 
17 Response latency (s) 
18 Distractor filtering index 
19x Loading weight ratios 
20 Delay 1 selectivity index 
21 Delay 2 selectivity index 
22x Delay 1 & 2 selectivity index (mean)
23x testing
%}

%% SELECT MEASURES TO BE CONSIDERED, REMOVE NEURONS NOT SHARED ACROSS MEASURE
nNeurons=length(funcMeasures{1});
allMeasures=[2 3 17 9 16 18 4 5 8 10 11 20 21]; %do not change this if you want to select specific measures. Do it in wantedMeasures
wantedMeasures=allMeasures; %[2 3 17 9 16 18 4 5 8 10 11 20 21]; %can be equal to allMeasures, or a subset, e.g. [2 3]
funcMeasures=keepSharedNeurons(funcMeasures,sort(allMeasures)); %this makes sure all measures have the same neurons, will cull some neurons!
%% Init
%% FIT MODELS AND IDENTIFY BEST FIT FOR EACH MEASURE
if HPC==0
    for nRegions=1:length(wantedRegions)
        nRegionStr=wantedRegions{nRegions};

        switch nRegionStr
            case {'3'}
                wantedMonkeys=2; %only one animal has 3 regions
        end

        for Monkey=wantedMonkeys
            MonkyID=['M' num2str(Monkey)];

            %Initialize breakpoint storage cell, pMC storage cell
            dDLPFCMeasuresCont={};
            vDLPFCMeasuresCont={};
            dDLPFCsurrPMeasuresCont={};
            vDLPFCsurrPMeasuresCont={};

            for measure=wantedMeasures
                %% Pre-process statistics

                % Measures which are proportions
                if measure<=8 %give each of responsive/selective neurons count of 1
                    wantedStat=funcMeasures{measure};

                    % convert index of nms neurons to a nNeuronsx1 binary matrix
                    desiredStatistic=zeros(nNeurons,1);
                    for neuron=nonan(setdiff(wantedStat,0))'
                        desiredStatistic(neuron)=1;
                    end

                % Measures which are not proportions
                elseif measure>=9
                    wantedStatIdx=funcMeasures{1};
                    wantedStat=funcMeasures{measure};
                    desiredStatistic=zeros(nNeurons,1);
                    i=1;
                    for neuron=wantedStatIdx'
                        desiredStatistic(neuron)=wantedStat(i);
                        i=i+1;
                    end
                end

                %% Fit gradients here
                [dDLPFCCont,vDLPFCCont,dDLPFCfitsPCell,vDLPFCfitsPCell]=fitGradients(folderPathEMPlot,folderPathNeuronIdx,electrodeMapping,...
                    desiredStatistic,measure,AP,Monkey,nRegionStr,loadFilteredData);

                for i=1:3
                    dDLPFCMeasuresCont{measure,i}=dDLPFCCont{i};
                    vDLPFCMeasuresCont{measure,i}=vDLPFCCont{i};
                    dDLPFCsurrPMeasuresCont{measure,i}=dDLPFCfitsPCell{i};
                    vDLPFCsurrPMeasuresCont{measure,i}=vDLPFCfitsPCell{i};
                end
            end

            %save breakpoints
            save([folderPathEM nRegionStr 'regions_' MonkyID 'breakpoints.mat'],'dDLPFCMeasuresCont','vDLPFCMeasuresCont')

            if Monkey==1
                dDLPFCsurrPMeasuresCont1=dDLPFCsurrPMeasuresCont;
                vDLPFCsurrPMeasuresCont1=vDLPFCsurrPMeasuresCont;
            elseif Monkey==2
                dDLPFCsurrPMeasuresCont2=dDLPFCsurrPMeasuresCont;
                vDLPFCsurrPMeasuresCont2=vDLPFCsurrPMeasuresCont;
            end
        end
        %% IDENTIFY BEST FITTING MODEL ACROSS MEASURES
        dDLPFCsurrPMeasuresCont2(3,:)=[]; %disselpct is invalid
        dDLPFCsurrPMeasuresContBoth=[dDLPFCsurrPMeasuresCont1;dDLPFCsurrPMeasuresCont2];
        vDLPFCsurrPMeasuresContBoth=[vDLPFCsurrPMeasuresCont1;vDLPFCsurrPMeasuresCont2];
        plotSurrPDist(dDLPFCsurrPMeasuresContBoth,vDLPFCsurrPMeasuresContBoth,3)
    end
    
elseif HPC==1
    for nRegions=1:length(wantedRegions)
        nRegionStr=wantedRegions{nRegions};

        switch nRegionStr
            case {'3'}
                wantedMonkeys=2; %only one animal has 3 regions
        end

        for Monkey=wantedMonkeys
            for measure=wantedMeasures
                %% Pre-process statistics
                % Measures which are proportions
                if measure<=8 %give each of responsive/selective neurons count of 1
                    wantedStat=funcMeasures{measure};

                    % convert index of nms neurons to a nNeuronsx1 binary matrix
                    desiredStatistic=zeros(nNeurons,1);
                    for neuron=nonan(setdiff(wantedStat,0))'
                        desiredStatistic(neuron)=1;
                    end

                % Measures which are not proportions
                elseif measure>=9
                    wantedStatIdx=funcMeasures{1};
                    wantedStat=funcMeasures{measure};
                    desiredStatistic=zeros(nNeurons,1);
                    i=1;
                    for neuron=wantedStatIdx'
                        desiredStatistic(neuron)=wantedStat(i);
                        i=i+1;
                    end
                end

                %% Fit gradients here
                fitGradientsHPC(folderPathEMPlot,folderPathNeuronIdx,electrodeMapping...
                ,desiredStatistic,measure,AP,Monkey,uniqueID);
            end
        end
    end
end