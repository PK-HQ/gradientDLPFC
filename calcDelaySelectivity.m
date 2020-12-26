function [neuronSelIndexes] =  calcDelaySelectivity(datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,assignLabel,delayPeriodBins)
nNeurons=size(datasets,1);

neuronSelIndexes=[];
for neuron = 1:nNeurons
    sessNumber = sessions(neuron,1);
    neuronData = squeeze(datasets(neuron,1:length(trials(sessNumber).val),:));
    targetLocations = AssignTrialLabel(trials(sessNumber).val,assignLabel);
    nLocations=length(unique(targetLocations));
    cont7Locations=[];
    for location = 1:nLocations
        locationTrials = find(targetLocations==location);
        locationData= neuronData(locationTrials,delayPeriodBins);
        meanLocationData=mean(locationData,2); %mean across trials
        meanLocationData=mean(meanLocationData(:)); %mean across bins
        cont7Locations=[cont7Locations,meanLocationData];
    end
    bestLocation=max(cont7Locations);
    worstLocation=min(cont7Locations);
    selectivityIndex=(bestLocation-worstLocation)/(bestLocation+worstLocation);
    neuronSelIndexes=[neuronSelIndexes;selectivityIndex];
end
    
        

end

