function [distractorFilteringRatioDT,distractorFilteringRatioTD] =  calcDistractorFiltering(datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,assignLabel)
targetPeriodBins=7:11;
distractorPeriodBins=33:37;
baselinePeriodBins=1:5;
nNeurons=size(datasets,1);

distractorFilteringRatioDT=[];
distractorFilteringRatioTD=[];
for neuron = 1:nNeurons
    sessNumber = sessions(neuron,1);
    neuronData = squeeze(datasets(neuron,1:length(trials(sessNumber).val),:));
    
    targetLocations = AssignTrialLabel(trials(sessNumber).val,assignLabel);
    nLocationsTar=length(unique(targetLocations));
    cont7LocationsTar=[];
    
    disLocations = AssignTrialLabel(trials(sessNumber).val,2);
    nLocationsDis=length(unique(disLocations));
    cont7LocationsDis=[];    
    
    
    
    for location = 1:nLocationsTar
        locationTrials = find(targetLocations==location);
        locationDataTar= neuronData(locationTrials,targetPeriodBins);
        meanLocationDataTar=mean(locationDataTar,2); %mean across trials
        meanLocationDataTar=mean(meanLocationDataTar(:)); %mean across bins
        cont7LocationsTar=[cont7LocationsTar,meanLocationDataTar];
    end
    [bestMeanTarFR,bestLoc]=max(cont7LocationsTar);

    for location = bestLoc
        locationTrials = find(disLocations==location);
        locationDataDis= neuronData(locationTrials,distractorPeriodBins);
        meanLocationDataDis=mean(locationDataDis,2); %mean across trials
        meanLocationDataDis=mean(meanLocationDataDis(:)); %mean across bins
    end
    
    distractorFilteringIdxDT=100*(meanLocationDataDis/bestMeanTarFR);
    distractorFilteringIdxTD=100*(bestMeanTarFR/meanLocationDataDis);

    distractorFilteringRatioDT=[distractorFilteringRatioDT;distractorFilteringIdxDT];
    distractorFilteringRatioTD=[distractorFilteringRatioTD;distractorFilteringIdxTD];
end
    
        

end

