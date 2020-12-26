function [relabelledArray]=removeGroundElectrodes(arrayXY)
%% Removes non-recording ground electrodes and relabels the recording electrodes with 1:32 
numElectrodes=size(arrayXY,1);

if numElectrodes==36
    %remove non-recording ground electrodes 1 9 28 and 36
    arrayXY(find(arrayXY(:,3)==1),:)=[];
    arrayXY(find(arrayXY(:,3)==9),:)=[];
    arrayXY(find(arrayXY(:,3)==28),:)=[];
    arrayXY(find(arrayXY(:,3)==36),:)=[];
    
    %relabel col3 electrodes with 1:32, as removal of ground electrodes
    %messes this up
    arrayXY(:,3)=numElectrodes-4:-1:1;
    relabelledArray=arrayXY;
elseif numElectrodes==18
    %remove non-recording ground electrodes 1 and 18
    arrayXY(find(arrayXY(:,3)==1),:)=[];
    arrayXY(find(arrayXY(:,3)==18),:)=[];
    
    %relabel col3 electrodes with 1:32, as removal of ground electrodes
    %messes this up
    arrayXY(:,3)=numElectrodes-2:-1:1;
    relabelledArray=arrayXY;
end
end
