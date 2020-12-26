function [folderPathCTDCond,fileSuffix]=getFilenameSuffix(testCondition,Trial_Label,folderPathCTD)
if strcmp(testCondition,'Error') || strcmp(testCondition,'error')
    if strcmp(Trial_Label,'Target') || strcmp(Trial_Label,'target')
        folderPathCTDCond=[folderPathCTD 'ET/'];
        fileSuffix='_ET';
    elseif strcmp(Trial_Label,'Distractor') || strcmp(Trial_Label,'distractor')
        folderPathCTDCond=[folderPathCTD 'ED/'];
        fileSuffix='_ED';
    elseif strcmp(Trial_Label,'Reward') || strcmp(Trial_Label,'reward')
        folderPathCTDCond=[folderPathCTD 'ER/'];
        fileSuffix='_ER';
    end
    
elseif strcmp(testCondition,'Correct') || strcmp(testCondition,'correct')
    if strcmp(Trial_Label,'Target') || strcmp(Trial_Label,'target')
        folderPathCTDCond=[folderPathCTD 'CT/'];
        fileSuffix='_CT';
    elseif strcmp(Trial_Label,'Distractor') || strcmp(Trial_Label,'distractor')
        folderPathCTDCond=[folderPathCTD 'CD/'];
        fileSuffix='_CD';
    elseif strcmp(Trial_Label,'Reward') || strcmp(Trial_Label,'reward')
        folderPathCTDCond=[folderPathCTD 'CR/'];
        fileSuffix='_CR';
    end
end
end

