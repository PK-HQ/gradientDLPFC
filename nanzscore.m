function [y,ymu,ysigma]=nanzscore(y)
if any(isnan(y(:)))
    ymu=nanmean(y);
    ysigma=nanstd(y);
    y=(y-repmat(ymu,length(y),1))./repmat(ysigma,length(y),1);
else
    [y,ymu,ysigma]=zscore(y);
end
end