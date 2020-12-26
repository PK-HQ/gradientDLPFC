function [loadingWeightRatio]=calcLoadingWeights(folderPathLW,folderPathNeuronIdx,folderPathRawData,maxbins,datasets,trials,bins_overlap,sessions,datasets_e,etrials,binsD1,binsD2)
lwRatioFilename=[folderPathLW 'lwRatio.mat'];
d1d2meanFilename=[folderPathLW 'd1d2mean.mat'];
decorrOutputFilename=[folderPathLW 'decorrOutputLW.mat'];
if exist(lwRatioFilename, 'file')==2
    load(lwRatioFilename)
else
    if exist(d1d2meanFilename, 'file')==2
        load(d1d2meanFilename)
    else
        [d1mean, d2mean]=getMeanD1D2(folderPathRawData,folderPathNeuronIdx,maxbins,...
            datasets,trials,bins_overlap,sessions,datasets_e,etrials,binsD1,binsD2);
        save(d1d2meanFilename,'d1mean', 'd2mean');
    end
    %NaN silent neurons
    d1silent=sum(d1mean==0,2)==7;
    d2silent=sum(d2mean==0,2)==7;
    d12silent=d1silent+d2silent;
    silentNeurons=find(d12silent>0);
    d1mean(silentNeurons,:)=NaN;
    d2mean(silentNeurons,:)=NaN;
    newNeuronIdx=1:632;
    if exist(decorrOutputFilename, 'file')==2
        load(decorrOutputFilename)
    else
        [Mcomp,Pcomp,a1,b1,minmi1,mi] = Decorrelation(d1mean,d2mean);
        save(decorrOutputFilename,'Mcomp','Pcomp','a1','b1','minmi1','mi');
    end

    %% GET LOADING WEIGHTS (Memory subspace/Motor Preparation subspace)
    N = size(d1mean,1); % number of neurons

    %construct the orthogonal basis for each subspace
    mem_space = gramschmidt(Mcomp,Mcomp(:,1));
    prp_space = gramschmidt(Pcomp,Pcomp(:,1));

    con_mem = zeros(N,1);
    con_prp = zeros(N,1);

    % contribution is the projection length of each neuron in the subspace
    for n = 1:N
        con_mem(n) = norm(mem_space(n,:),2);
        con_prp(n) = norm(prp_space(n,:),2);
    end

    % get loading weight ratio 
    % i.e., contribution to Memory subspace/Motor Preparation subspace
    loadingWeightRatio=con_mem./con_prp;

    %save
    save(lwRatioFilename,'loadingWeightRatio');
end

end