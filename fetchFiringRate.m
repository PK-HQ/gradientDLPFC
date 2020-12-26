function [firingrateLocCorrect,firingrateLocError]=fetchFiringRate(...
    trials,etrials,neuronSession,neuronSessionIdx,locationCoord,firingrateSessionCorrect,firingrateSessionError)

%get correct and error trial indices

sessionDataCorrect=[trials(neuronSession).val];
sessionDataError=[etrials(neuronSession).val];

correctLoc=AssignTrialLabel(sessionDataCorrect,1);
idxCorrectLoc=find(correctLoc==locationCoord);

errorLoc=AssignTrialLabel(sessionDataError,1);
idxerrorLoc=find(errorLoc==locationCoord);

%get firing rates for trial index of desired location
firingrateLocCorrect=firingrateSessionCorrect(1,idxCorrectLoc,:);
firingrateLocError=firingrateSessionError(1,idxerrorLoc,:);


end