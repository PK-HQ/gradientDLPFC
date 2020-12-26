function removeDuplicatesCombineMonkeyCells(folderPathNeuronIdx,versionStr,anatLocations)
M1Neu=load([folderPathNeuronIdx,'M1_regionalNeuronsIdx_' versionStr '.mat']);
M1Neu=struct2mat(M1Neu);
M2Neu=load([folderPathNeuronIdx 'M2_regionalNeuronsIdx_' versionStr '.mat']);
M2Neu=struct2mat(M2Neu);
load([folderPathNeuronIdx 'duplications.mat']);
idx=find(dup==1);
nReg=size(M1Neu,2);
regionalNeuronsIdx=zeros(632,nReg);

for regCol=1:nReg
    %pool neurons
    M1RegNeu=nonzeros(M1Neu(:,regCol));
    M2RegNeu=nonzeros(M2Neu(:,regCol));
    combRegNeu=[M1RegNeu;M2RegNeu];
    
    %rm duplicate
    [~,pos]=intersect(combRegNeu,idx);
    combRegNeu(pos)=[];
    
    %put into zero padded mat
    nRegNeu=numel(combRegNeu);
    fprintf('%s: %.0f\n',anatLocations{regCol},nRegNeu)
    regionalNeuronsIdx(1:nRegNeu,regCol)=sort(combRegNeu);
end
    
numel(nonzeros(regionalNeuronsIdx))

save(['data/neuronidx/regionalNeuronsIdx' versionStr '.mat'],'regionalNeuronsIdx')
end