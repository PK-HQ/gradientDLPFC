function [datasetsCorrect, datasetsError, etrials, trials, countCorrectIdx, countErrorIdx] = extractCorrectErrorTrials(datasets, sessions, fulltrials)
datasetFull=datasets;
nSessions=unique(sessions)';
countCorrectIdx=0;
countErrorIdx=0;

%get first and last trials number of each session
reference=[];
for session=nSessions
    reference(session,1)=min(find(sessions==session));
    reference(session,2)=max(find(sessions==session));
end

datasetsCorrect=[];
datasetsError=[];

% create correct and error trials structure (session names)
for session=nSessions
    trials(session).name=fulltrials(session).name;
    etrials(session).name=fulltrials(session).name;
end

for session=nSessions
    %extract neuron indexes grouped by session, e.g. neurons for session 1
    sessionIdx=sessions==session;
    %extract neurons for session, e.g. session 1
    sessionData=datasetFull(sessionIdx,:,:);
    %extract trial indexes grouped by correct/error, e.g. correct and wrong trials
    trialCorrectOrError={fulltrials(session).val.failure};
    cond=cellfun(@isempty,trialCorrectOrError);
    find(cond);
    correctIdx=find(cond);
    errorIdx=find(~cond);
    %extract correct and wrong trials
    allTrials=1:size(trialCorrectOrError,2);
    trialIdxCorrect=allTrials(correctIdx);
    trialIdxError=allTrials(errorIdx);
    %create correct and wrong trial trials struct
    trials(session).val=fulltrials(session).val(correctIdx);
    etrials(session).val=fulltrials(session).val(errorIdx);
    %recompile correct trial datasets
    ncols1=1:size(trialIdxCorrect,2);
    rowStart=reference(session,1);
    rowEnd=reference(session,2);
    datasetsCorrect(rowStart:rowEnd,ncols1,:)=sessionData(:,trialIdxCorrect,:);
    %recompile error trial datasets
    ncols2=1:size(trialIdxError,2);
    rowStart=reference(session,1);
    rowEnd=reference(session,2);
    datasetsError(rowStart:rowEnd,ncols2,:)=sessionData(:,trialIdxError,:);
    
    %% etrials
    trialsStruct=fulltrials(session).val;
    %errorROWS=trialsStruct(errorIdx);
    %etrials(session).val=errorROWS;
    countCorrectIdx=countCorrectIdx+numel(correctIdx);
    countErrorIdx=countErrorIdx+numel(errorIdx);
end
end








