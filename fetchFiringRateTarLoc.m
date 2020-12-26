function [firingrateLocT]=fetchFiringRateTarLoc(...
        neuron,trials,etrials,session,locations,locationCoord,firingrateSessionCorrect,firingrateSessionError,correctTrials,correctAndErrorTrials,minTrials)
if correctAndErrorTrials==0 & correctTrials==1 || correctAndErrorTrials==1 & correctTrials==1 || correctAndErrorTrials==1 & correctTrials==0
    sessionData=[trials(session).val];

    %label CT locs
    TLoc=AssignTrialLabel(sessionData,1);
    idxTLoc=find(TLoc==locationCoord);

    %get firing rates for trial index of desired location
    firingrateLocT=firingrateSessionCorrect(1,idxTLoc,:);

elseif correctAndErrorTrials==0 & correctTrials==0
    sessionData=[etrials(session).val];

    %label CT locs
    TLoc=AssignTrialLabel(sessionData,1);
    idxTLoc=find(TLoc==locationCoord);

    %get firing rates for trial index of desired location
    firingrateLocT=firingrateSessionError(1,idxTLoc,:);

    
end


end