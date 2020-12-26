function [d1mean, d2mean]=getMeanD1D2(folderPath,folderPathNeuronIdx,maxbins,datasets,trials,bins_overlap,sessions,datasets_e,etrials,winLP11,winLP22)
%% make a 632 x 3 matrix of neuron 632 idx | session | session's neuron idx
sessNeurons=[];
for i=unique(sessions)'
    nNeurons=sum(sessions==i);
    tmp=1:nNeurons;
    sessNeurons=[sessNeurons,tmp];
end
sessNeurons=sessNeurons';
neurons632=1:632;
majorIdx=[neurons632',sessions,sessNeurons];

for neuron=majorIdx(:,1)' %wantedNeurons
    fprintf('%.0f\n',neuron)
    neuronSession=majorIdx(neuron,2);
    neuronSessionIdx=majorIdx(neuron,3);

    %get firing rates for this neuron
    correctLocsFR=datasets(neuron,:,:);
    errorLocsFR=datasets_e(neuron,:,:);

    %get firing rates per location, put into a 1x7 cell
    for location=1:7
        
        %get trials for wanted location 
        [correctLocFR,errorLocFR]=fetchFiringRate(trials,etrials,neuronSession,neuronSessionIdx,location,correctLocsFR,errorLocsFR);
        
        baselineFR=(correctLocFR(:,:,1:6));
        baselineMeans=nanmean(baselineFR,3);%baselineMean=nanmean(baselineMeans,2);
        
        d1FR=correctLocFR(:,:,winLP11);
        d1Means=nanmean(d1FR,3)-baselineMeans;d1Mean=nanmean(d1Means,2);
        d1mean(neuron,location)=d1Mean;
        
        d2FR=correctLocFR(:,:,winLP22);
        d2Means=nanmean(d2FR,3)-baselineMeans;d2Mean=nanmean(d2Means,2);
        d2mean(neuron,location)=d2Mean;
    end
end
end























