function [binsTar,binsDis,binsD1,binsD2]=defEventTimebins(bins_overlap)
%gives indices of bins for each event of Target and Distractor
%presentation, delays 1 and 2
binCenters=mean(bins_overlap,1);
targetDuration=300; D1Duration=1000; distractorDuration=300; D2Duration=1000;
midD1=targetDuration+D1Duration/2;
endD1=targetDuration+D1Duration;
midD2=targetDuration+D1Duration+distractorDuration+D2Duration/2;
endD2=targetDuration+D1Duration+distractorDuration+D2Duration;

startTar=0;
endTar=targetDuration;
startDis=targetDuration+D1Duration;
endDis=targetDuration+D1Duration+distractorDuration;
binsD1=find(binCenters==midD1,1):find(binCenters==endD1,1); %last 500ms delay 1
binsD2=find(binCenters==midD2,1,'first'):find(binCenters==endD2,1,'first'); %last 500ms delay 1
binsTar=find(bins_overlap(1,:)==startTar,1,'first'):find(bins_overlap(2,:)==endTar,1,'first'); %300ms target presentation
binsDis=find(bins_overlap(1,:)==startDis,1,'first'):find(bins_overlap(2,:)==endDis,1,'first'); %300ms distractor presentation
end