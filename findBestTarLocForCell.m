function [maxOfBin]=findBestTarLocForCell(...
    firingRateLocCT)
%firingRateMeanCorrectAll(1:30,1)=NaN;
%firingRateMeanErrorAll(1:30,1)=NaN;
%firingRateSpreadCorrectAll(1:30,1)=NaN;
%firingRateSpreadErrorAll(1:30,1)=NaN;
%for timeWin=1:maxbins
%get FR for a 100ms timebin
dataC=firingRateLocCT(1,:,7:11);
maxOfBin=nanmean(dataC(:));
%store this mean and spread data into a 1x1xbins container that contains all 5 timebins



end