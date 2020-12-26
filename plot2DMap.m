function plot2DMap(folderPathEMPlot,folderPathNeuronIdx)
%% Load single neuron data
load([folderPathNeuronIdx 'regionalNeuronsIdxBothMonkeys.mat']);
load([folderPathNeuronIdx 'electrodeMappingCurved.mat']);

%% Actual plotting
allNeurons=sort(nonzeros(reshape(regionalNeuronsIdx,numel(regionalNeuronsIdx),1))); %1:632 array
plotElectrodeMap(folderPathEMPlot,electrodeMapping,allNeurons);

end
