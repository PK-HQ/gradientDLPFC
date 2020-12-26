function bootstrapAICCont=fitDisModelIterHPC_shuffle(x,yActual,allBreaks,yLab,xLab,statLimits,statSpacing,...
    statName,wantedregion,regionColor,regionStr,folderName,uniqueID,Monky)
%per monkey id
switch Monky
    case {1}
        MonkyInitial='P';
    case {2}
        MonkyInitial='J';
end

%parameters
maxiters=100;
n2seg=0;
n3seg=0;
statOfInterest='surrogateR2';
robust=1;
bootstrapAICCont=nan(maxiters,1);
bootstrapR2Cont{maxiters,5}=[];
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
if strcmp(regionStr,'iLPFC')
    breakpoints=[b2];
end
regionBoundaries=[b1 b2;
                  b3+0.01 b4;
                  b5+0.01 b6];
jitterBoundary=1.5;

xFull=x;
yFull=yActual;
%% Get linear model fit
%linMdl=fitlm(xFull(:),yFull(:),'RobustOpts',robustStr);
rng('shuffle','twister')
rng('shuffle','twister')
for iter=1:maxiters
    %% Generate shuffle data 
    x=xFull(:);    
    y=yFull(randperm(length(yFull)));
    bootstrapR2Cont{iter,4}=x;
    bootstrapR2Cont{iter,5}=y;
    %% Find best 1 and 2 breakpoints
    %calculate each break from the actual data, and residual of the resultant 2 lines flanking it
    best2BreakPoints=find2Breaks(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    best1BreakPoint=find1Break(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    
    %% 1 segment model (surrogate)
    mdlA=fitlm(x,y,'RobustOpts',robustStr);
    switch statOfInterest
        case {'AIC'}
            %AIC
            SSR_1Seg = sum((mdlA.Residuals.Raw).^2);
            K_1Seg=mdlA.NumEstimatedCoefficients;
            N_1Seg=mdlA.NumObservations;
            AICc1=AIC(SSR_1Seg,N_1Seg,K_1Seg);
            bootstrapAICCont(iter,1)=AICc1;
        case {'RMSE'}
            [arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);
            bootstrapR2Cont{iter,1}=arsq1whole;
            %RMSE1Seg=sqrt(mdl.SSE/mdl.DFE);
        case {'surrogateR2'}
            [arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);
            bootstrapR2Cont{iter,1}=arsq1whole;
    end
    
    %% 2 segment model (shuffle actual y values within the region's bounds)
    xb1 =best1BreakPoint(1,1);
    if ~isnan(xb1)==0
        validBreak=0;
    else
        n2seg=n2seg+1;
        edges = [min(x),xb1,max(x)];
        bins = discretize(x,edges);
        coeffs=[];
        [idx1,~]=find(bins == 1);
        [idx2,~]=find(bins == 2);
        validBreak=(numel(idx1)>2)+(numel(idx2)>2)+sum(isnan(edges)); %check if empty and no-nan
    end
    switch validBreak
        case {2}
            %m1
            x1=x(idx1);
            y1=y(idx1);
            %rand shuffle y
            %y1 = y1(randperm(length(y1)));
            mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
            coeffs=[coeffs;mdlA.Coefficients.Estimate(1);mdlA.Coefficients.Estimate(2)];
            %m2
            x2=x(idx2);
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
                    bootstrapAICCont(iter,2)=AICc2;
                case {'RMSE'}
                    [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
                    %RMSE2Seg=sqrt((mdlA.SSE + mdlB.SSE)/(mdlA.DFE + mdlB.DFE));
                case {'surrogateR2'}
                    [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
                    bootstrapR2Cont{iter,2}=arsq2whole;
            end
        otherwise
            AICc2=NaN;
            arsq2whole=NaN;
    end

    %% 3 segment model (shuffle actual y values within the region's bounds)
    switch regionStr
        case {'midV'}
            %don't fit
            arsq3whole=NaN;
            AICc3=NaN;
        case {'midD'}
            %fit model
            xb1 =min(best2BreakPoints);
            xb2 =max(best2BreakPoints);
            if ~isnan(xb1)==0 || ~isnan(xb1)==0
                validBreak=0;
            else
                n3seg=n3seg+1;
                edges = [min(x),xb1,xb2,max(x)];
                bins = discretize(x,edges);
                coeffs=[];
                [idx1,~]=find(bins == 1);
                [idx2,~]=find(bins == 2);
                [idx3,~]=find(bins == 3);
                validBreak=(numel(idx1)>2)+(numel(idx2)>2)+(numel(idx3)>2)+sum(isnan(edges));
                %validBreak=~isempty(idx1)+~isempty(idx2)+~isempty(idx3)+sum(isnan(edges));
            end
            switch validBreak
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
                            bootstrapAICCont(iter,3)=AICc3;
                        case {'RMSE'}
                            [arsq3whole,~,~,~]=getMdlStats(mdlA,mdlB,mdlC);
                             %RMSE3Seg=sqrt((mdlA.SSE + mdlB.SSE + mdlC.SSE)/(mdlA.DFE + mdlB.DFE + mdlC.DFE));
                        case {'surrogateR2'}
                            [arsq3whole,~,~,~]=getMdlStats(mdlA,mdlB,mdlC);
                            bootstrapR2Cont{iter,3}=arsq3whole;
                    end
                otherwise
                    AICc3=NaN;
                    arsq3whole=NaN;
            end
    end
    if n2seg>maxiters & n3seg>maxiters
        break
    end
    %}
end
save(['/home/svu/lsitpk/HPC_PK/data/em/gradient/shuffled/shuffle_' MonkyInitial '_' regionStr '_' folderName '_' num2str(maxiters) 'i_' uniqueID '.mat'],'bootstrapR2Cont');
    
end
