
function [nms,lms,f_nms_interact,f_nms_mainEffect,d1_sel,d2_sel,tar_sel,dis_sel,f_del1,f_del2,size_RF_d2,size_RF_d1,size_RF_dT,cs] = TwoWayAnova(datasets,trials,sessions,bins_overlap)
winBaseline=1:6;
winLPT=7:11;
winLPDis=33:37;
winLP11=22:32;
winLP22=48:58;
%datapath='/Users/PK/Google Drive/MI6/L-lab/software/matlab/pfc/script/Code-for-Parthasarathy-et-al-master/datasets_aishu.mat';
%load(datapath);
bins_overlap=bins_overlap(:,1:59);
datasets=datasets(:,:,1:60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs ANOVA on neural data in spike counts defined by the
% input variable bins_overlap to identify NMS, LMS and CS neurons as defined in the
% paper. Further, you can also compute the receptive field size for the
% neuron's responsiveness during the trial, f-stats for NMS and delay 1 and
% delay 2 selectivity. The neurons identified by this code and the
% F-stats computed here are used in Fig 5 and Fig 6 in Parthasarathy et al.
% Any questions?? Please contact Aishwarya Parthasarathy at aishu.parth@gmail.com
% 30th August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input -
% datasets: Size - Nneurons x Ntrials x Nbins.
% This matrix contains the activity of all the neurons used in decoding
% during correct trials in relevant time bins_overlap. Ntrials is the maximum
% number of trials recorded by the neurons. If a neuron is recorded with
% less than Ntrials trials, then the rest of the values in the matrix till
% Ntrials is zero-padded. Nbins is the length of the time bins_overlap for which
% this decoding was performed. Therefore the 5th neuron's activity in the
% 6th trial during the 8th time bin is denoted by datasets(5,6,8)
%
% trials:  Size - array of structs - 1 x Nsession
% Nsession is the number of recorded sessions from which we pool neurons to
% form datasets. Nsession also equals the number of unique values in the
% first column of the sessions array. Each struct is in turn an array of
% structs containing all the timing and stimulus information presented in
% all the correct trials in that particular sessions.
%
% sessions: Size - Nneurons x 2
% Each element in the first column refers to the sessions in which the
% neuron was recorded. Also, this value as the index of the trials variable
% fetches all the trials performed in the sessions where the said neuron was
% recorded.
%
% bins_overlap: Size - 2 X Nbins
% The first row specifies the starting time point of the bin and the second
% row specifies the ending time point.
%
% Output -
% nms - Nonlinear Mixed Selective -  Indices of the nms neurons amongst the
% neurons described in datasets.
%
% lms - Linear Mixed Selective -Indices of the lms neurons amongst the
% neurons described in datasets.
%
% cs - Classically Selective - Indices of the cs neurons amongst the neurons
% described in datasets
%
% f_nms - f-stats for the 2-way interaction term for each of the nms neurons.
% Same length as nms.
%
% d1_sel - indices of neurons that exhibit target selectivity during Delay
% 1 (750 - 1250 ms)
%
% d2_sel - indices of neurons that exhibit target selectivity during Delay
% 2 (2050 - 2550 ms)
%
% f_del1 - f-stats from ANOVA for those neurons identified in d1_sel.
%
% f_del2 - f_stats from ANOVA for those neurons identified in d2_sel.
%
% size_RF_d2 - no of responsive locations for each neuron during Delay 2
% (2050 - 2550ms)
%
% size_RF_d1 - no of responsive locations for each neuron during Delay 1
% (750 - 1250 ms).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the variables
h = zeros(size(datasets,1),3);
f_stat_nms_interact = zeros(size(datasets,1),1);
f_stat_nms_mainEffect = zeros(size(datasets,1),1);
f_stat_del1 = zeros(size(datasets,1),1);
f_stat_del2 = zeros(size(datasets,1),1);
size_RF_dT = zeros(size(datasets,1),1);
size_RF_d1 = zeros(size(datasets,1),1);
size_RF_d2 = zeros(size(datasets,1),1);
p1 = zeros(size(datasets,1),1);
p2 = zeros(size(datasets,1),1);
p3 = zeros(size(datasets,1),1);
% Looping through all the neurons in datasets
for i = 1:size(datasets,1)
    % tr1 is the sessions in which neuron i was recorded
    tr1 = sessions(i,1);
    % tmp contains all the neural data that was recorded during correct
    % trials in sessions tr1 and is stored as spike counts in time bins_overlap
    % defined by the variable bins_overlap. n_trials denoted the number of correct
    % trials in sessions tr1.
    n_trials = length(trials(tr1).val);
    tmp = squeeze(datasets(i,1:n_trials,:));
    %tar period
    dT = mean(tmp(:,winLPT),2);
    dDis = mean(tmp(:,winLPDis),2);
    
    % d1 is the mean neural activity for delay 1 during all correct trials recorded for
    % neuron i
    %d1 = mean(tmp(:,find(bins_overlap(2,:)==750):find(bins_overlap(2,:)==1250)),2);
    d1 = mean(tmp(:,winLP11),2);
    % d2 is the mean neural activity for delay 2 during all correct trials
    % recorded for neuron i
    %d2 = mean(tmp(:,find(bins_overlap(2,:)==1950):find(bins_overlap(2,:)==2450)),2);
    d2 = mean(tmp(:,winLP22),2);
    % base defines the mean baseline activity defined by the activity 300
    % ms prior to the target onset
    base = mean(tmp(:,winBaseline),2);
    % Avg delay activity across all trials.
    avg_d1 = mean(d1);
    avg_d2 = mean(d2);
    % Subtracting the mean baseline activity from the delay 1 and delay 2
    % activity.
    
    % SHOULD NOT SUBTRACT!!!
    %d1 = (d1 - mean(base));
    %d2 = (d2 - mean(base));
    
    % If the average activity in delay 1 or delay 2 is 0, the test for
    % selectivity is skipped. This is to avoid low firing neurons.
    if avg_d1==0 || avg_d2==0
        continue
    end
    % tar is a vector containing all the target labels of the correct trials
    % performed by the animal in the sessions that the neuron was recorded
    tar = AssignTrialLabel(trials(tr1).val,1);
    dis = AssignTrialLabel(trials(tr1).val,2);
    
    % Consolidating delay 1 and delay 2 to feed to anovan (pre-processing
    % anovan)
    dat = [d1;d2];
    % g1 is a vector that identifies whether said neural activity comes
    % from d1 or d2. 1000 if its d1 and 2000 if its d2. Its just a random
    % placeholder to define delay 1 and delay 2 identity (pre-processing).
    g1 = [repmat(1000,length(tar),1);repmat(2000,length(tar),1)];
    % g2 is a vector for target labels for all the values in dat (pre-processing)
    g2 = [tar';tar'];
    % Performing 2 way anova with g1 and g2 as factors. h stores the 3
    % p-values for each neuron, one for each factor and one for the two way
    % interaction.
    [h(i,:),tabD12] = anovan(dat,{g1 g2},'model','full','display','off');

    % f-stats for the two-way main effect
    f_stat_nms_mainEffect(i,1) = tabD12{[3],6};

    % f-stats for the two-way interaction
    f_stat_nms_interact(i,1) = tabD12{4,6};
    % Performing anova to check for target selectivity in Tar
    [p3(i),table3,~] = anova1(dT,tar,'off');
    % Performing anova to check for target selectivity in Dis
    [p4(i),table4,~] = anova1(dDis,dis,'off');
    % Performing anova to check for target selectivity in Delay 1
    [p1(i),table1,~] = anova1(d1,tar,'off');
    % f-stats for target selectivity in Delay 1
    f_stat_del1(i,1) = table1{2,5};
    % Initiliazing the size of Receptive Field (target) to be 0
    % Looping through all target locations
    for iu = 1:length(unique(tar))
        % If the neuron is responsive to a target location, increase the
        % counter on the size of the Receptive Field.
        ind_tar = tar==iu;
        % A ttest is performed comparing the activity for a particular
        % target location to the baseline to checck for responsiveness
        [~,p_ttest] = ttest(dT(ind_tar),mean(base));
        if p_ttest < 0.05
            size_RF_dT(i) = size_RF_dT(i) + 1;
        end
    end
    % Initiliazing the size of Receptive Field (target) to be 0
    % Looping through all target locations
    for iu = 1:length(unique(tar)) %for each target
        % If the neuron is responsive to a target location, increase the
        % counter on the size of the Receptive Field.
        ind_tar = tar==iu; %get trials for that target
        % A ttest is performed comparing the activity for a particular
        % target location to the baseline to checck for responsiveness
        %[~,p_ttest] = ttest(d1(ind_tar),mean(base));
        [~,p_ttest] = ttest(d1(ind_tar),mean(base)); %compare d1 activity for target vs baseline
        if p_ttest < 0.05 %if SD, add one to counter
            size_RF_d1(i) = size_RF_d1(i) + 1;
        end
    end
    % Similar check for target selectivity in Delay 2
    [p2(i),table2,~] = anova1(d2,tar,'off');
    % f-stats for target selectivity in Delay 2
    f_stat_del2(i,2) = table2{2,5};
    % Initializing and computing the size of receptive field for delay 2
    for iu = 1:length(unique(tar))
        ind_tar = tar==iu;
        [~,p_ttest] = ttest(d2(ind_tar),mean(base));
        if p_ttest < 0.05
            size_RF_d2(i) = size_RF_d2(i) + 1;
        end
    end
end


% h(:,1) - p-value for the factor defining the effect due to task epoch
% h(:,2) - p-value for the factor defining the effect of target locations. 
% h(:,3) - p-value for the interaction term

cs = find(h(:,3)>0.05 & h(:,1)>0.05 & h(:,2) < 0.05); % Identifies Classical selectivity
lms = find(h(:,3) > 0.05 & h(:,1) < 0.05 & h(:,2) < 0.05 ); %Identifies LMS
nms = find(h(:,3) < 0.05 & h(:,3) ~= 0); % Identifies NMS
tar_sel = find(p3 < 0.05); % Neurons that are selective in Tar
dis_sel = find(p4 < 0.05); % Neurons that are selective in Dis
d1_sel = find(p1 < 0.05); % Neurons that are selective in Delay 1
d2_sel = find(p2 < 0.05); % Neurons that are selective in Delay 2
%f_nms_mainEffect = f_stat_nms_mainEffect(nms,:); 
%f_nms_interact = f_stat_nms_interact(nms); % F-stats for NMS neurons (Fig 5f)
f_nms_mainEffect = f_stat_nms_mainEffect; 
f_nms_interact = f_stat_nms_interact; % F-stats for NMS neurons (Fig 5f)
%f_del1 = f_stat_del1(d1_sel); %F-stats for delay 1 and delay 2 selectivity
%f_del2 = f_stat_del2(d2_sel);
f_del1 = f_stat_del1; %F-stats for delay 1 and delay 2 selectivity
f_del2 = f_stat_del2;
fprintf('============================\n\nSelectivity type: NMS=%d, LMS=%d, CS=%d, Total=%d\n', numel(nms), numel(lms), numel(cs), numel(nms)+numel(lms)+numel(cs))
fprintf('Selective neurons: D1=%d, D2=%d, D1 & D2=%d\n\n============================\n\n', numel(d1_sel), numel(d2_sel), numel(d1_sel)+numel(d2_sel)-numel(intersect(d1_sel,d2_sel)))
%sprintf('Selectivity f-stats: NMS=%d\nD1=%d\nD2=%d', f_nms,f_del1,f_del2')
%sprintf('RF: D1=%d\nD2=%d', size_RF_d1',size_RF_d2')
end







