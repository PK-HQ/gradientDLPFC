function [cleanedData,outlierIdx]=rmoutlier(data,mode)
switch mode
    case {'mean'} %mean and 3 s.d.
        stdDev = nanstd(data(:)); % Compute standard deviation
        meanValue = nanmean(data(:)); % Compute mean
        outlierIdx = abs(data-meanValue) > (3* stdDev);% Create a binary map of outlier position
        data(outlierIdx)=[];
        cleanedData=data;
    case {'median'} %median and 3 median abs dev
        outlierIdx=abs(data - median(data)) > (3*mad(data,1));
        data(outlierIdx)=[];
        cleanedData=data;
end
end