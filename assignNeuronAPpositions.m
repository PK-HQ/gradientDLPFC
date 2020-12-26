function electrodeAPpositions=assignNeuronAPpositions(electrodeMapping, neuronElectrodeIDStrs)
%% find neuron's electrode in an array with all electrode x- and y- positions, assign the neuron its AP and DV positions
%neuronElectrodeIDStrs=cell array w neuron IDs and their array and electrode string, i.e. {5, g12, c2}
%electrodeMapping= electrode string, and electrode x- and y-coordinates after rotation and translation, 
%i.e. {g12, 1, c2, 2.41, 3.01} 

%set params
electrodeAPpositions={};
electrodeIDCol=1;
electrodeArrayStrCol=2;
electrodeIDStrCol=3;

%for each neuron assign it an x- y-coord based on its electrode ID str
for n=1:length(neuronElectrodeIDStrs)
    electrodeID=neuronElectrodeIDStrs{eachNeuron,electrodeIDCol}; %get this electrode's array str
    electrodeArrayStr=neuronElectrodeIDStrs{eachNeuron,electrodeArrayStrCol}; %get this electrode's array str
    electrodeIDStr=neuronElectrodeIDStrs{eachNeuron,electrodeIDStrCol}; %get this electrode's str
    
    wantedElectrode=find(electrodeMapping(:,1)==electrodeArrayStr & electrodeMapping(:,4)==electrodeIDStr); %find electrode str in array
    
    [x,y]=electrodeMapping{wantedElectrode,2:3}; %get electrode's x- and y-
    
    electrodeAPpositions{eachNeuron,1}=electrodeID; %save in new array
    electrodeAPpositions{eachNeuron,2}=x;
    electrodeAPpositions{eachNeuron,3}=y;
end
end