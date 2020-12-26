function calcFuncMeasures(folderPathNeuronIdx,folderPathLW,folderPathRawData,...
    datasets,trials,bins_overlap,sessions,datasets_e,etrials,anatLocations,...
    trialLabel,binsD1,binsD2,maxbins)
%% CALCULATES ALL MEASURES USED IN THE GRADIENT ANALYSIS (+ EXTRAS)

%% Define params and containers
%define containers
nNeurons=size(datasets,1);
stimulusAppearanceBins=find(bins_overlap(1,:)>=0,1,'first'):size(bins_overlap,2);
neuronCount=1:nNeurons';
testing=ones(nNeurons,1);

%% CALCULATE BULK STATS cs Lms nms rfs tarsel d1sel d2sel f_d1 f_d2 nms-interact nms-main
[nms,Lms,f_nms_interact,f_nms_mainEffect,...
    d1_sel,d2_sel,tar_sel,dis_sel,f_del1,f_del2,size_RF_d2,size_RF_d1,size_RF_dT,cs]...
    = TwoWayAnova(datasets,trials,sessions,bins_overlap);
%{
[nmsE,LmsE,f_nms_interactE,f_nms_mainEffectE,...
    d1_selE,d2_selE,tar_selE,dis_selE,f_del1E,f_del2E,size_RF_d2E,size_RF_d1E,size_RF_dTE,csE]...
    = TwoWayAnova(datasets_e,etrials,sessions,bins_overlap);
%}
    
%% CALCULATE SELECTIVITY INDEX (T/D1/D2)
[selIdxD2] =  calcDelaySelectivity(...
    datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,trialLabel,binsD1);
[selIdxD1] =  calcDelaySelectivity(...
    datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,trialLabel,binsD2);
[selIdx] =  tarSelectivity(datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,trialLabel);
meanSelIdxD1D2=nanmean([selIdxD1,selIdxD2],2);

%% CALCULATE DISTRACTOR FILTERING
[distractorFilteringRatioDT,distractorFilteringRatioTD] =  calcDistractorFiltering(...
    datasets,trials,sessions,bins_overlap,nNeurons,stimulusAppearanceBins,trialLabel);
%}
%% CALCULATE LATENCY (MS)
offset=-2;
baselineBins=1:find(bins_overlap(1,:)>=0,1,'first')-1;
trialBinTotalTar=find(bins_overlap(1,:)>=0,1,'first'):size(bins_overlap,2)+offset;

[respNeuronTotal,respNeuronD1,respNeuronD2,earliestRespBin,tarLatency]=getRespNeuronAndLatency(...
    nNeurons,sessions,datasets,trials,datasets_e,etrials,anatLocations,...
    baselineBins,trialBinTotalTar,binsD1,binsD2,1,bins_overlap);

%% CALCULATE LOADING WEIGHTS MEMORY/MOTOR PREP
lwRatio=calcLoadingWeights(folderPathLW,folderPathNeuronIdx,folderPathRawData,maxbins,datasets,...
    trials,bins_overlap,sessions,datasets_e,etrials,binsD1,binsD2);

%% COMPILE FUNCTIONAL MEASURES AND SAVE, (COL 1) WITH LABELS (COL 2)
funcMeasures={neuronCount;tar_sel;dis_sel;d1_sel;d2_sel;cs;Lms;nms;...
    size_RF_dT;size_RF_d1;size_RF_d2;f_nms_interact;f_nms_mainEffect;...
    f_del1;f_del2(:,2);selIdx;tarLatency;distractorFilteringRatioDT;lwRatio;...
    selIdxD1;selIdxD2;meanSelIdxD1D2;testing};
for i=1:length(funcMeasures)
    if size(funcMeasures{i},1)<size(funcMeasures{i},2)
        funcMeasures{i}=funcMeasures{i}';
    end
    if length(funcMeasures{i})<nNeurons
        funcMeasures{i}=make632(funcMeasures{i});
    elseif length(funcMeasures{i})>nNeurons
        fprintf('Error, incorrect measure size, too many neurons')
    end
end

%save
funcMeasures(:,2)={'Neuron index (1:632)';'Target selective';'distractor selective';...
    'Delay 1 selective';'Delay 2 selective';'Classical selective';...
    'Linear mixed selective';'Nonlinear mixed selective (NMS)';'Receptive field size';...
    'D1 field size';'D2 field size';'NMS interaction term';'NMS main effect';...
    'Delay 1 selective f-stat';'Delay 2 selective f-stat';'Selectivity index';...
    'Response latency (s)';'Distractor filtering index';'Loading weight ratios';...
    'Delay 1 selectivity index';'Delay 2 selectivity index';...
    'Delay 1 & 2 selectivity index (mean)';'testing'};
save([folderPathNeuronIdx 'functionalMeasures.mat'],funcMeasures);
end