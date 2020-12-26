function [best2BreakPoints,tiedBreakpoints]=find2Breaks(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr)
RMSECell={};
bestBreakPointCont=[];

for bb=1:numel(breakpoints)
    breakpointMin=breakpoints(bb)-jitterBoundary;
    breakpointMax=breakpoints(bb)+jitterBoundary;
    if breakpointMax>max(x)
        breakpointMax=max(x);
    end
    stepSize=0.1;
    RMSECont=[];
    for breakpoint=breakpointMin:stepSize:breakpointMax

        %get data for first 2 regions
        startRegion1=regionBoundaries(bb,1);
        endRegion2=regionBoundaries(bb+1,2);
        dataIdx = x >= startRegion1 & x <= endRegion2;
        xData=x(dataIdx);
        yData=y(dataIdx);
        
        %get data before break
        dataIdx1 = xData <= breakpoint;
        xData1=xData(dataIdx1);
        yData1=yData(dataIdx1);
        
        %get data after break
        dataIdx2 = xData > breakpoint;
        xData2=xData(dataIdx2);
        yData2=yData(dataIdx2);
        
        if numel(xData1)<3 || numel(xData2)<3 %numel(xData1)==0 || numel(xData2)==0
            RMSEwhole=NaN;
        else
            %fit linear, get coefficients
            mdl1=fitlm(xData1,yData1,'RobustOpts',robustStr);
            mdl2=fitlm(xData2,yData2,'RobustOpts',robustStr);
            RMSEwhole=sqrt((mdl1.SSE + mdl2.SSE)/(mdl1.DFE + mdl2.DFE));

        end
        RMSECont=[RMSECont; breakpoint RMSEwhole];
    end
    RMSECell{bb,1}=RMSECont;
end

% Get 2 best breaks for region boundaries
bbb=1;
for bb=1:numel(breakpoints)
    breakpointResiduals=RMSECell{bb,1};
    if isempty(breakpointResiduals)
        minResidual=NaN;
        minRMSEBreakPoints=NaN;
        breakpt=NaN;
    else
        minResidual=min(breakpointResiduals(:,2));
        breakpt=breakpoints(bb);
        if numel(minResidual)>1
            %remove duplicate minimum residual
            minResidual=unique(minResidual);
        end
        minRMSEBreakPoints=breakpointResiduals(breakpointResiduals(:,2)==minResidual,1);
    end
    %take closest break to anat boundary
    [~,idx]=min(abs(minRMSEBreakPoints-breakpt));
    bestBreakPoint=minRMSEBreakPoints(idx);
    bestBreakPointCont(bb,1:2)=[minResidual,bestBreakPoint];
    tiedBreakpoints{bb}=minRMSEBreakPoints;
end

bestBreakPointCont(:,3)=2;
if sum(~isnan(bestBreakPointCont(:,1)))==0
    best2BreakPoints=NaN;
else
    bestOf2Breaks=bestBreakPointCont(bestBreakPointCont(:,1)==min(bestBreakPointCont(:,1)),2);
    bestBreakPointCont(bestBreakPointCont(:,2)==bestOf2Breaks,3)=1;
    best2BreakPoints=bestBreakPointCont(:,2);
end
end