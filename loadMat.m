function [dataMat]=loadMat(folderPath,fileName)
dataStruct=load([folderPath fileName]);
dataMat=cell2mat(struct2cell(dataStruct));
end