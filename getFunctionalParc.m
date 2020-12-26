function getFunctionalParc(folderPathNeuronIdx,anatLocations)
% load electrode map
load([folderPathNeuronIdx 'electrodeMappingCurved.mat']);
load([folderPathNeuronIdx 'regionalNeuronsIdx.mat']);
allNeurons=sort(nonzeros(reshape(regionalNeuronsIdx,numel(regionalNeuronsIdx),1)));
getRegionalIdx(folderPathNeuronIdx,electrodeMapping,allNeurons,anatLocations);
end