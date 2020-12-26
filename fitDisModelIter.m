function AdjR2Cont=fitDisModelIter(x,yActual,allBreaks,yLab,xLab,statLimits,statSpacing,statName,...
    wantedregion,regionColor,regionStr,folderName,maxiters)
%maxiters=1000;
statOfInterest='RMSE';
robust=0;
AdjR2Cont=nan(maxiters,3);
if robust==1
    robustStr='on';
    robustPlotStr='_R';
elseif robust==0
    robustStr='off';
    robustPlotStr=[];
end

sigTestStatStr='Spearman';

%conts
b1=allBreaks(1);
b2=allBreaks(2);
b3=allBreaks(3);
b4=allBreaks(4);
b5=allBreaks(5);
b6=allBreaks(6);
breakpoints=[b2 b4];
regionBoundaries=[b1 b2;
                  b3+0.01 b4;
                  b5+0.01 b6];
jitterBoundary=1.5;

%% Get stats only
for iter=1:maxiters
    tic;
    yShuffled2region=yActual(randperm(length(yActual)));

    if strcmp(regionStr,'sLPFC')
        twoRegionElectrodes=find(x<=regionBoundaries(5));
        x2region=x(twoRegionElectrodes);
        yShuffled2region=yShuffled2region(twoRegionElectrodes);
    elseif strcmp(regionStr,'iLPFC')
        x2region=x;
        yShuffled2region=yShuffled2region;
    end

    %calculate each break from the actual data, and residual of the resultant 2 lines flanking it
    best2BreakPoints=find2Breaks(x,yShuffled,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    best1BreakPoint=find1Break(x2region,yShuffled2region,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    %% 1 segment model (shuffle all y values)
    %fit model
    x2region = x(:); % easier as column vectors
    y = yShuffled2region(:);
    
    %rand shuffle y
    %y = y(randperm(length(y)));
    mdlA=fitlm(x2region,y,'RobustOpts',robustStr);
    switch statOfInterest
        case {'AIC'}
            %AIC
            SSR_1Seg = sum((mdlA.Residuals.Raw).^2);
            K_1Seg=mdlA.NumEstimatedCoefficients;
            N_1Seg=mdlA.NumObservations;
            AICc1=AIC(SSR_1Seg,N_1Seg,K_1Seg);
        case {'RMSE'}
            [arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);
            %RMSE1Seg=sqrt(mdl.SSE/mdl.DFE);
    end
    
    
    %% 2 segment model (shuffle actual y values within the region's bounds)
    %fit model
    xb1 =best1BreakPoint(1,1);
    edges = [min(x2region),xb1,max(x2region)];
    bins = discretize(x2region,edges);
    coeffs=[];
    [idx1,~]=find(bins == 1);
    [idx2,~]=find(bins == 2);

    
    nEmpty=~isempty(idx1)+~isempty(idx2);
    switch nEmpty
        case {2}
            %m1
            x1=x2region(idx1);
            y1=y(idx1);
            %rand shuffle y
            %y1 = y1(randperm(length(y1)));
            mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
            coeffs=[coeffs;mdlA.Coefficients.Estimate(1);mdlA.Coefficients.Estimate(2)];
            %m2
            x2=x2region(idx2);
            y2=y(idx2);
            %rand shuffle y
            %y2 = y2(randperm(length(y2)));
            mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
            coeffs=[coeffs;mdlB.Coefficients.Estimate(1);mdlB.Coefficients.Estimate(2)];
            %AIC
            switch statOfInterest
                case {'AIC'}
                    SSR_2Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2);
                    K_2Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients;
                    N_2Seg=mdlA.NumObservations+mdlB.NumObservations;
                    AICc2=AIC(SSR_2Seg,N_2Seg,K_2Seg);
                case {'RMSE'}
                    [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
                    %RMSE2Seg=sqrt((mdlA.SSE + mdlB.SSE)/(mdlA.DFE + mdlB.DFE));
            end
        otherwise
            arsq2whole=NaN;
    end

    %% 3 segment model (shuffle actual y values within the region's bounds)
    %fit model
    x = x(:); % easier as column vectors
    y = yShuffled(:);
    xb1 =min(best2BreakPoints);
    xb2 =max(best2BreakPoints);
    edges = [min(x),xb1,xb2,max(x)];
    bins = discretize(x,edges);
    coeffs=[];
    [idx1,~]=find(bins == 1);
    [idx2,~]=find(bins == 2);
    [idx3,~]=find(bins == 3);
    nEmpty=~isempty(idx1)+~isempty(idx2)+~isempty(idx3);
    switch nEmpty
        case {3}
            %m1
            x1=x(idx1);
            y1=y(idx1);
            %rand shuffle y
            %y1 = y1(randperm(length(y1)));
            mdlA=fitlm(x1,y1,'RobustOpts',robustStr);

            %m2
            x2=x(idx2);
            y2=y(idx2);
            %rand shuffle y
            %y2 = y2(randperm(length(y2)));
            mdlB=fitlm(x2,y2,'RobustOpts',robustStr);

            %m3
            x3=x(idx3);
            y3=y(idx3);
            %rand shuffle y
            %y3 = y3(randperm(length(y3)));
            mdlC=fitlm(x3,y3,'RobustOpts',robustStr);
            switch statOfInterest
                case {'AIC'}
                    %AIC
                    SSR_3Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2) + sum((mdlC.Residuals.Raw).^2);
                    K_3Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients+mdlC.NumEstimatedCoefficients;
                    N_3Seg=mdlA.NumObservations+mdlB.NumObservations+mdlC.NumObservations;
                    AICc3=AIC(SSR_3Seg,N_3Seg,K_3Seg);
                case {'RMSE'}
                    [arsq3whole,~,~,~]=getMdlStats(mdlA,mdlB,mdlC);
                     %RMSE3Seg=sqrt((mdlA.SSE + mdlB.SSE + mdlC.SSE)/(mdlA.DFE + mdlB.DFE + mdlC.DFE));
            end
        otherwise
            arsq3whole=NaN;
    end
    %% Define reference AIC
    switch statOfInterest
        case {'AIC'}
            refAICc=min([AICc1,AICc2,AICc3]);
        case {'RMSE'}
            AdjR2Cont(iter,1:numel([arsq1whole,arsq2whole,arsq3whole]))=[arsq1whole,arsq2whole,arsq3whole];
    end
    toc;
    fprintf('%.0f (%.1f s)\n',iter,toc);
end
save(['/Volumes/Users/PK/Desktop/HPC_PK/data/em/gradient/shuffled/' regionStr '_' folderName '100.mat'],'AdjR2Cont');
    
end
