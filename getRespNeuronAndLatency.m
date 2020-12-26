function [respNeuronTotal,respNeuronD1,respNeuronD2,earliestRespBin,tarLatency]=getRespNeuronAndLatency(...
    nNeurons,sessions,datasets,trials,...
    datasets_e,etrials,anatLocations,...
    baselineBins,trialBinTotal,trialBinD1,trialBinD2,assignLabel,bins_overlap)
% get best location for each neuron (highest firing rate)
bestLocationForCellCont=getBestLocation(folderPathNeuronIdx,datasets,trials,...
            sessions,datasets_e,etrials,anatLocations,0,0,1);

respNeuronTotal=[];
respNeuronD1=[];
respNeuronD2=[];
earliestRespBin_7loc=nan(632,8);
tarLatency=nan(632,1);

for neuron = 1:nNeurons %for each neuron
    tr1 = sessions(neuron,1);
    neuronTrialAllLocs = squeeze(datasets(neuron,1:length(trials(tr1).val),:));
    baseline=mean(mean(neuronTrialAllLocs(:,baselineBins)));
    targetLabelForEachTrial = AssignTrialLabel(trials(tr1).val,assignLabel);
    for targetLocation = 1:length(unique(targetLabelForEachTrial)) %for each target location
        %get trial indexes of a specific target location
        trialsOfSpecificTargetLocation = find(targetLabelForEachTrial==targetLocation);
        
        %total responsive
        for timebin=trialBinTotal
            %timebin
            %type 1: consecutive bins are sig
            %{
            bin1=neuronTrialAllLocs(trialsOfSpecificTargetLocation,timebin);
            bin2=neuronTrialAllLocs(trialsOfSpecificTargetLocation,timebin+1);
            [~,bin1pSig] = ttest(bin1, baseline);
            [~,bin2pSig] = ttest(bin2, baseline);
            a = bin1pSig < 0.05;
            b = bin2pSig < 0.05;
            if a==1 && b==1
                respNeuron=[respNeuron,i];
            end
            %}
            %type 2: joined bins are sig

            bin1=neuronTrialAllLocs(trialsOfSpecificTargetLocation,timebin);
            [~,binpSig] = ttest(bin1, baseline);
            a = binpSig < 0.05;
            if a==1
                respNeuronTotal=[respNeuronTotal,neuron];
                time=mean(bins_overlap(:,timebin));
                earliestRespBin_7loc(neuron,targetLocation,timebin)=time;
                %if D1
                if timebin>=trialBinD1(1) && timebin<=trialBinD1(2)
                    respNeuronD1=[respNeuronD1,neuron];
                %if D2
                elseif timebin>=trialBinD2(1) && timebin<=trialBinD2(2)
                    respNeuronD2=[respNeuronD2,neuron];
                end
            end

        end

    end
end

%% Use the neuron's preferred stimulus location if possible, else pick the earliest resp timebin
earliestRespBin_7loc(earliestRespBin_7loc==0)=NaN;
for i=1:numel(bestLocationForCellCont)
    if ~isnan(bestLocationForCellCont(i))
        neuronFR=earliestRespBin_7loc(i,bestLocationForCellCont(i),:);
        earliestBin=min(neuronFR,[],3);
    else
        neuronFR=earliestRespBin_7loc(i,:,:);
        earliestBin=min(neuronFR,[],3);
        earliestBin=min(earliestBin,[],2);
    end
    tarLatency(i)=earliestBin;
end

%convert 0 to nan
earliestRespBin_7loc(earliestRespBin_7loc==0)=NaN;
%get min across timebins (across pages)
earliestRespBin=min(earliestRespBin_7loc,[],3);
%get min across lcoations (across columns)
earliestRespBin=min(earliestRespBin,[],2);


respNeuronTotal=unique(respNeuronTotal);
respNeuronD1=unique(respNeuronD1);
respNeuronD2=unique(respNeuronD2);
end