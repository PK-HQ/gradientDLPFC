function [best1BreakPoint,minRMSEBreakPoints]=find1Break(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr)
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
        startRegion1=min(regionBoundaries(:));
        endRegion3=max(regionBoundaries(:));
        xData=x;
        yData=y;
        
        %get data before break
        dataIdx1 = x >= startRegion1 & x <= breakpoint;
        xData1=xData(dataIdx1);
        yData1=yData(dataIdx1);
        
        %get data after break
        dataIdx2 = x > breakpoint & x <= endRegion3;
        xData2=xData(dataIdx2);
        yData2=yData(dataIdx2);
        
        if numel(xData1)<3 || numel(xData2)<3
            RMSE=NaN;
            RMSEwhole=NaN;
        else
            %fit linear, get coefficients
            mdl1=fitlm(xData1,yData1,'RobustOpts',robustStr);
            mdl2=fitlm(xData2,yData2,'RobustOpts',robustStr);

            %store RMSE to optimize
            %RMSE1=sum((mdl1.RMSE).^2);
            %RMSE2=sum((mdl2.RMSE).^2);
            %RMSE=RMSE1+RMSE2;
            %[RMSEwhole,~,~,~]=getMdlStats(mdl1,mdl2,[]);
            RMSEwhole=sqrt((mdl1.SSE + mdl2.SSE)/(mdl1.DFE + mdl2.DFE));

        end
        RMSECont=[RMSECont; breakpoint RMSEwhole];

    end
    RMSECell{bb,1}=RMSECont;
end

% Get 1 best breaks for region boundaries
bbb=1;
for bb=1:numel(breakpoints)
    breakpointResiduals=RMSECell{bb,1};
    
    if isempty(breakpointResiduals)
        minResidual=NaN;
        minRMSEBreakPoints=NaN;
        breakpt=NaN;
    else
        minResidual=min(breakpointResiduals(:,2));
        breakpt=breakpoints(bbb);
        if numel(minResidual)>1
            %take the one nearest to 0
            minResidual=unique(minResidual);
        end
        minRMSEBreakPoints=breakpointResiduals(breakpointResiduals(:,2)==minResidual,1);
    end
    
    

    if numel(minResidual)>1
        %take the one nearest to 0
        minResidual=unique(minResidual);
    end
    %take closest break to anat boundary
    [~,idx]=min(abs(minRMSEBreakPoints-breakpt));
    bestBreakPoint=minRMSEBreakPoints(idx);
    bestBreakPointCont(bb,1:2)=[minResidual,bestBreakPoint];
    bbb=bbb+1;
    tiedBreakpoints{bb}=minRMSEBreakPoints;
end

%% Rank breakpoints, tag with 1 or 2
bestBreakPointCont(:,3)=2;

if sum(~isnan(bestBreakPointCont(:,1)))==0
    best1BreakPoint=NaN;
else
    bestOf2Breaks=bestBreakPointCont(bestBreakPointCont(:,1)==min(bestBreakPointCont(:,1)),2);
    best1BreakPoint=bestOf2Breaks(1);
end


end