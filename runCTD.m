function performance = runCTD(datasets,m,st,sessions,trialLabel,trials,...
    bins_overlap,etrials,datasets_e,m_e,st_e,Test,parcellatedNeuronsIdx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function decodes a trial label from the datasets built from the
% recorded neurons. Specifically, the function evaluates the performance of
% the decoder when you dissociate the  time windows from which you pick the
% training and testing data. This datasets may be built from a pseudo-population or
% from simultaneously recorded ensemble of neurons. This code generates
% results for Fig 2a, 2e, 4a, 4f, 6a, 6b, 6c, 6e, 6f and 6g. This
% code is run 1000 times to generate a distribution of performances. The
% mean of this distribution of performances is plotted as a heatmap
% in these figures.
% Any questions?? Please contact Aishwarya Parthasarathy at aishu.parth@gmail.com
% 30th August 2017
%%%%%%%%%%%%%%%%%

% Debugging guide:
% Are your data structures correct?
% did you feed correct matrices for correct and error datasets and not the combined dataset?
% are parcellation matrices the correct n-parcels and size?
% is bin size correct?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs-

% datasets : Size - Nneurons x Ntrials x Nbins.
% This matrix contains the activity of all the neurons used in decoding
% during correct trials in relevant time bins_overlap. Ntrials is the maximum
% number of trials recorded by the neurons. If a neuron is recorded with
% less than Ntrials trials, then the rest of the values in the matrix till
% Ntrials is zero-padded. Nbins is the length of the time bins_overlap for which
% this decoding was performed. Therefore the 5th neuron's activity in the
% 6th trial during the 8th time bin is denoted by datasets(5,6,8)

% m : Size - Nneurons x 1
% Each element in this array is the mean activity across all trials during
% the baseline period (300ms prior to the target onset) of a neuron

% st: Size - Nneurons x 1
% Each element in this array is the standard deviation of the activity
% activity across all trials during the baseline period (300ms prior to the
% target onset) of a neuron

% sessions: Size - Nneurons x 1
% Each element in the first column refers to the sessions in which the
% neuron was recorded. Also, this value as the index of the trials variable
% fetches all the trials performed in the sessions where the said neuron was
% recorded.

% trialLabel - 'target' or 'distractor'
% this string refers to the trial label to be decoded from the datasets.

% trials: Size - struct 1 x Nsessions
% Nsessions is the number of recorded sessionss from which we pool neurons to
% form datasets. Nsessions also equals the number of unique values in the
% first column of the sessionss array.

% dataset_e,m_e,st_e,etrials are the same as datasets,m,st and trials
% repectively but for error trials.

% Test is a string that takes the value 'Correct' or 'Error' to denote if
% we are decoding correct or error trials. This string is useful in
% differentiating the decoding in Fig 2 and Fig 4 of Parthasarathy et al
% Note while decoding the target from data recorded during error trials,
% the decoder is trained using data from correct trials and tested using data
% from error trials.
%
% Output -
% performance: Size - Nbinstrain x Nbinstest
% An element (x,y) in this array refers to the performance of the decoder when
% trained using data from the xbin and tested using data from ybin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('/Users/PK/Desktop/software/matlab/pfc/dataset/datasets_tar.mat');

%% Builds a set of training and testing trial labels to build the training and testing datasets
% that will eventually be fed to the decoder.
[trainingLabel,trainingTrials,testingLabel,testingTrials] = MakeTrialSet(trialLabel,trials,etrials,sessions,Test);

%% Initializing the training and testing datasets with zeros
nNeurons=size(datasets,1);
nTrainingTrials=size(trainingTrials,2);
nTestingTrials=size(testingTrials,2);
nBins=size(datasets,3);
trainingData = zeros(nNeurons,nTrainingTrials,nBins);
testingData = zeros(nNeurons,nTestingTrials,nBins);

%diagnostics
%trainingDataStats=zeros(1:size(datasets,1),3);
%testingDataStats=zeros(1:size(datasets,1),3);

%% Filling up training data with datasets data (nNeuron X 1500 Trials X Bins)
% Filling up train_data from datasets. All the values are z_scored using m,st to normalize the firing rate across neurons
countEmptyRowsTraining=0;
for neuron = 1:size(datasets,1)
    neuronTrialFR=(datasets(neuron,trainingTrials(neuron,:),:));
    neuronBaselineMean=m(neuron);
    neuronBaselineSTD=st(neuron);
    zscoredData=(neuronTrialFR-neuronBaselineMean)./neuronBaselineSTD;
    if sum(isnan(zscoredData(:)))>0 || sum(isinf(zscoredData(:)))>0 %remove zscore=NaN (i.e. 0/0) or =inf (i.e. 10/0)
        %fprintf('NAN | Neuron: %d Mean:%d STD:%d \n', neuron, neuronBaselineMean, neuronBaselineSTD)
        countEmptyRowsTraining=countEmptyRowsTraining+1;
        continue %dont fill this row, skip to next neuron
    end
    trainingData(neuron,:,:) = zscoredData; %zscore correct dataset data
end
%size(trainingData) %see some rows of nans

%% Filling up testing data with datasets data (nNeuron X 100 Trials X Bins)
% Filling up test_data from datasets. All the values are z_scored using m,st to normalize the firing rate across neurons
%fprintf('\ntest_data %d x %d x %d\n',size(test_data))
%fprintf('\ntest_trials %d x %d\n',size(test_trials))
%fprintf('\ntest_trials %d\n',size(datasets,1))
countEmptyRowsTesting=0;
if strcmp(Test,'Correct') || strcmp(Test,'correct')
    %fprintf('Correct trials (testing set)')
    for neuron = 1:size(datasets,1)
        neuronTrialFR=datasets(neuron,testingTrials(neuron,:),:);
        neuronBaselineMean=m(neuron);
        neuronBaselineSTD=st(neuron);
        testingData(neuron,:,:) = (neuronTrialFR-neuronBaselineMean)./neuronBaselineSTD; %zscore correct dataset data
    end
elseif strcmp(Test,'Error') || strcmp(Test,'error')
    %fprintf('Error trials (testing set)')
    for neuron = 1:size(datasets,1)
        neuronTrialFR=datasets_e(neuron,testingTrials(neuron,:),:);
        neuronBaselineMean=m_e(neuron);
        neuronBaselineSTD=st_e(neuron);
        zscoredData=(neuronTrialFR-neuronBaselineMean)./neuronBaselineSTD; %zscore correct dataset data
        if sum(isnan(zscoredData(:)))>0 || sum(isinf(zscoredData(:)))>0 %remove zscore=NaN (i.e. 0/0) or =inf (i.e. 10/0)
            %fprintf('NAN | Neuron: %d Mean:%d STD:%d \n', neuron, neuronBaselineMean, neuronBaselineSTD)
            countEmptyRowsTesting=countEmptyRowsTesting+1;
            continue %dont fill this row, skip to next neuron
        end
        testingData(neuron,:,:) = zscoredData; 
        %testingData(neuron,:,:) = (neuronTrialFR-neuronBaselineMean)./neuronBaselineSTD; %zscore error dataset data
    end
    
else
    %fprintf('MISSING TEST STRING!!!')
    for neuron = 1:size(datasets,1)
        neuronTrialFR=datasets(neuron,testingTrials(neuron,:),:);
        neuronBaselineMean=m(neuron);
        neuronBaselineSTD=st(neuron);
        testingData(neuron,:,:) = (neuronTrialFR-neuronBaselineMean)./neuronBaselineSTD; %zscore correct dataset data
    end
end

%% Extract neurons from desired lpfc parcel
trainingData=trainingData(parcellatedNeuronsIdx,:,:); %extract only neurons for the desired lpfc parcel
testingData=testingData(parcellatedNeuronsIdx,:,:);
%figure()

%% For each training data timebin, denoise then decode testing data timebin
% Looping through the training time bins_overlap
for trainingBin = 1:size(bins_overlap,2)
    % Looping through the testing time bins_overlap
    for testingBin = 1:size(bins_overlap,2)
        % De-noising the datasets to feed into the decoder
        [trainingDataNew,testingDataNew] = Build_DataSet(squeeze(trainingData(:,:,trainingBin)),squeeze(testingData(:,:,testingBin)));
        % Decoding and computing the percentage of the trials that the
        % decoder classified correctly for corresponding training and
        % testing window
        [performance(trainingBin,testingBin) ] = ComputePerformance(trainingDataNew,testingDataNew,trainingLabel,testingLabel);
    end
end

end

%% Function to create a matrix of training and testing labels
function [trainingLabel,trainingTrials,testingLabel,testingTrials] = MakeTrialSet(trialLabel,trials,etrials,sessions,Test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds a matrix of trial indices and labels to build the
% training set. Trials from every sessions are split into training and
% testing pool. Further, two uniform distributions of 7 target labels was
% built to define the training and testing set's trial label. The length of
% these uniform distribution is defined by the number of trials used to
% train and test the decoder. For example, if the first trial label in the training set
% is target location 1, the function picks one trial for each neuron in the
% training pool with the target presented at target 1. Similarly, this is
% repeated for all the trial labels in the training and testing set.
% Inputs -
% trialLabel - string - 'target' or 'distractor'
% trials - struct- passed on from the main function
% etrials - struct - passed on from the main function
% sessions - Nneurons x 2 matrix - passed on from the main function.
% Outputs -
% trainingLabel - Ntraintrials x 1 matrix - uniform distribution of the 7
% trial_labels (target or distractor)
% trainingTrials - Nneurons x Ntraintrials - indices of trials with
% trial label specified by trainingLabel for all Nneurons.
% testingLabel and testingTrials are similar to trainingLabel and
% trainingTrials respectively but the trials to build testingTrials are
% chosen from the testing pool. Usually the length of trainingLabel is set to
% 1500 and the length of testingLabel is set to 100.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deciphering the trial label to decode and assigning label the
% corresponding value, label = 1 when target is decoded and label = 2 when
% distractor is decoded.
if strcmp(trialLabel,'Target') || strcmp(trialLabel,'target')
    label=1;
elseif strcmp(trialLabel,'Distractor') || strcmp(trialLabel,'distractor')
    label=2;
elseif strcmp(trialLabel,'reward') || strcmp(trialLabel,'reward')
    label=5;
end
% Dividing trials from each sessions into training and testing groups
for nSession = 1:length(trials)
    % Randomly picking 50% of the trials in a sessions to be under the
    % training pool
    train_num = randsample(length(trials(nSession).val),round(0.50*(length(trials(nSession).val))));
    % Assigning the other 50% to be the testing pool. While decoding error
    % trials it will be replaced with 1:length(etrials(i_sessions).val)
    if strcmp(Test,'Correct') || strcmp(Test,'correct')
        test_num = setdiff(1:length(trials(nSession).val),train_num);
    elseif strcmp(Test,'Error') || strcmp(Test,'error')
        test_num = 1:length(etrials(nSession).val);
    else
        test_num = setdiff(1:length(trials(nSession).val),train_num);
    end
    % Building the train_set with trial details for each sessions using
    % train_num
    train_set(nSession).val = trials(nSession).val(train_num);
    % Storing their original indices from trials
    train_set(nSession).orgind = train_num;
    % Storing the target labels of all the trials in train_set. Please note
    % that AssignTrialLabel function is specific for our datasets. This
    % function identifies the label (target or distractor location) for each 
    % trial. This needs to be modified if you are not using the datasets
    % used in Parthasarathy et al and subsequently the lines of code using
    % the output of AssignTrialLabel.
    train_set(nSession).tarlabel = AssignTrialLabel(train_set(nSession).val,label);
    % Similar variables for test_set
    if strcmp(Test,'Correct') || strcmp(Test,'correct')
        test_set(nSession).val = trials(nSession).val(test_num);
        test_set(nSession).orgind = test_num;
        test_set(nSession).tarlabel = AssignTrialLabel(test_set(nSession).val,label);
    elseif strcmp(Test,'Error') || strcmp(Test,'error')
        test_set(nSession).val = etrials(nSession).val;
        test_set(nSession).orgind = test_num;
        test_set(nSession).tarlabel = AssignTrialLabel(test_set(nSession).val,label);
    else
        test_set(nSession).val = trials(nSession).val(test_num);
        test_set(nSession).orgind = test_num;
        test_set(nSession).tarlabel = AssignTrialLabel(test_set(nSession).val,label);
    end
    train_num=[];test_num=[];
end
% Setting Ntraintrials
nTrainingTrials = 1500;
% Setting Ntesttrials
test_tr = 100;
% nNeuronInSession stores the number of neurons recorded in each sessions
% Initializing nNeuronInSession
nNeuronInSession = zeros(length(trials),1);
for nSession = 1:length(trials)
    nNeuronInSession(nSession,1) = length(find(sessions(:,1)==nSession));
end
% Creates a uniform distribution of trial labels (based on trialLabel)
% between 1 and 7 of length nTrainingTrials
if strcmp(Test,'Correct') || strcmp(Test,'correct')
    trainingLabel = randsample(1:7,nTrainingTrials,true);
    %Creates a uniform distribution of trial labels between 1 and 7 of length
    %test_tr
    testingLabel = randsample(1:7,test_tr,true);
elseif strcmp(Test,'Error') || strcmp(Test,'error')
    trainingLabel = randsample([2 3 5 6],nTrainingTrials,true);
    testingLabel = randsample([2 3 5 6],test_tr,true);
else
    trainingLabel = randsample([2 3 5 6],nTrainingTrials,true);
    testingLabel = randsample([2 3 5 6],test_tr,true);
end
% id is a counter for the number of cells used in this analysis.
nNeuron = 1;
% i_sessions loops through the number of recorded sessionss used in this
% analysis.
for nSession = 1:length(trials)
    % Checking if there are any neurons recorded in the sessions
    if nNeuronInSession(nSession,1)~=0
        % if there are neurons in that recorded sessions, loop through the
        % length of training set. For each value in trainingLabel, find the
        % trials from the training pool for that recorded sessions
        % (represented as i_sessions) And repeat this for
        % every trial label in trainingLabel
        for nTrainingSetVal = 1:nTrainingTrials %1:1500
            % Initializing temporary variables
            train_label_tmp=[];ind=[];
            % Finding all the trials in the training pool with the nth value of trainingLabel.
            ind = find(train_set(nSession).tarlabel==trainingLabel(nTrainingSetVal));
            % Sample with replacement from ind as many times as the number
            % of neurons in the recorded sessions (i_sessions)
            train_label_tmp = (randsample(length(ind),nNeuronInSession(nSession),true));
            % Build trainingTrials with the selected trials with
            % train_label_tmp. Note id is the cell counter and
            % id+nNeuronInSession(i_sessions,1) is the counter after adding all
            % the neurons recorded in i_sessions.
            a=train_set(nSession).orgind(ind(train_label_tmp));
            trainingTrials(nNeuron:nNeuron+nNeuronInSession(nSession,1)-1,nTrainingSetVal) = a;
        end
        % Go through the same loop for test trials to build hte testingTrials
        for nTrainingSetVal = 1:length(testingLabel)
            ind=[];test_label_tmp=[];
            ind = find(test_set(nSession).tarlabel==testingLabel(nTrainingSetVal));
            test_label_tmp = (randsample(length(ind),nNeuronInSession(nSession),true));
            testingTrials(nNeuron:nNeuron+nNeuronInSession(nSession)-1,nTrainingSetVal) = test_set(nSession).orgind(ind(test_label_tmp));
        end
    end
    nNeuron = nNeuron+nNeuronInSession(nSession);
end
clearvars -except trainingLabel trainingTrials testingLabel testingTrials
end

%% Function to build the training and testing set
% This function denoises the train and test datasets using PCA. The data
% used to decode is the data projected onto the principal components.
function [trainingDataNew,testingDataNew] = Build_DataSet(train_data,test_data)
% Initializing a matrix to store all the neuron indices with NaN values.
% These are neurons with very low firing.
ind_nan=[];
% Checking for neurons with NaN values within the train data
for i = 1:size(train_data,1)
    if ~isfinite(train_data(i,1))
        ind_nan = [ind_nan i];
    end
end
% Checking for neurons with NaN values within the test data
for i = 1:size(test_data,1)
    if ~isfinite(test_data(i,1))
        ind_nan = [ind_nan i];
    end
end
% Pick all the neurons that had NaNs in the train and/or test data.
ind_nan = unique(ind_nan);
% Getting rid of these neurons in the train and test datasets.
train_data(ind_nan,:)=[];test_data(ind_nan,:)=[];
% Creating a PCA space with training and testing data.
A = [train_data';test_data'] - mean([train_data';test_data'],2);
[V,D] = eig(cov(A));
%[coeff,score,latent] = princomp([train_data';test_data']);
score = A*fliplr(V);
% Computing the proportion of explained variance for each component
latent = sort(diag(D),'descend');
latent = cumsum(latent);
latent = latent/latent(end);
% Finding the number of components explaining 90% of the variance.
expl_var = dsearchn(latent,0.9);
% The denoised train data is the projection of the original data on the
% first n components of the PCA space created in line 215.
trainingDataNew = score(1:size(train_data,2),1:expl_var);
% Denoising the test data similarly.
testingDataNew = score(size(train_data,2)+1:end,1:expl_var);

end

%% Calculating decoding performance using an LDA
function [performance] = ComputePerformance(train_data,test_data,training_label,testing_label)
% class is the predicted target label from the test_data
[class,~,~,~,~] = classify(test_data,train_data,training_label);
% Checking how different the predicted label is from the actual label for
% all the test trials
perf = class-testing_label';
% Extracting all the correctly predicted target label
perf_ind = find(perf==0);
% The performance of the decoder is computed as the percentege of number of correct
% predictions in the decoding.
performance = length(find(perf==0))*100/length(testing_label);


end


