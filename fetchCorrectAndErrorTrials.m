function [datasets,trials,m,st,bins_overlap,sessions,regional,datasets_e,etrials,m_e,st_e]=...
    fetchCorrectAndErrorTrials(folderPathRawData,fileName,timeBins,saveFiles)
% separates correct and error trials (and corresponding data) from binned data.

%% LOAD RAW
fprintf('\n==================== Extract correct and error trials ====================\n')
fprintf('Loading raw data...')
datapath=[folderPathRawData fileName];
[datasets, bins_overlap, regional, sessions, fulltrials, m, st] = loadData(datapath, timeBins);
fprintf('Raw data loaded!\n')

%% SEPARATE CORRECT AND ERROR TRIALS
fprintf('Extracting correct and error trials...')
[datasetsCorrect, datasetsError, etrials, trials, countCorrectIdx, countErrorIdx] = extractCorrectErrorTrials(datasets, sessions, fulltrials);
datasets=datasetsCorrect;
datasets_e=datasetsError;
m_e=m;
st_e=st;
fprintf('Trials extracted! (%d correct and %d error trials)\n',countCorrectIdx,countErrorIdx)

%% SAVE FILES
switch saveFiles
    case {1}
        fprintf('Saving correct and error trials...')
        %normal naming set
        fileNameCorr=[folderPathRawData 'datasetsCorrectTrials.mat'];
        fileNameErr=[folderPathRawData 'datasetsErrorTrials.mat'];

        %Set A (for HPC)
        fileNameCorr_A=[folderPathRawData 'datasetsCorrectTrials1.mat'];
        fileNameErr_A=[folderPathRawData 'datasetsErrorTrials1.mat'];

        %Set B (for HPC)
        fileNameCorr_B=[folderPathRawData 'datasetsCorrectTrials2.mat'];
        fileNameErr_B=[folderPathRawData 'datasetsErrorTrials2.mat'];

        %Saving sets
        save(fileNameCorr,'datasets','trials','bins_overlap','regional','sessions','m','st');
        save(fileNameErr,'datasets_e','etrials','bins_overlap','regional','sessions','m_e','st_e');

        save(fileNameCorr_A,'datasets','trials','bins_overlap','regional','sessions','m','st');
        save(fileNameErr_A,'datasets_e','etrials','bins_overlap','regional','sessions','m_e','st_e');

        save(fileNameCorr_B,'datasets','trials','bins_overlap','regional','sessions','m','st');
        save(fileNameErr_B,'datasets_e','etrials','bins_overlap','regional','sessions','m_e','st_e');
        
        fprintf('Saved correct and error trials!\n')
end
end

