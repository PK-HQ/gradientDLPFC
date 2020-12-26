function [refAICc,pvalues,breakpointCellCont,arsqCont]=fitDiscontinuousModel(x,y,allBreaks,yLab,xLab,statLimits,statSpacing,statName,wantedregion,...
regionColor,shuffledAdjR2Cont,folderName,regionStr,wantedregion2,subplotIdx,Monky,robust)

%monkey id
switch Monky
    case {1}
        MonkyStr={''};
        MonkyInitial='P';
        shuffledAdjR2Cont(:,3)=nan(length(shuffledAdjR2Cont),1); %padding
    case {2}
        MonkyStr={''};
        MonkyInitial='J';
        if strcmp(regionStr,'vDLPFC')
            shuffledAdjR2Cont(:,3)=nan(length(shuffledAdjR2Cont),1); %padding
        end
end

breakpointCellCont=cell(1,3);
dotSize=100; %80;
AICcont=[];
statOfInterest='RMSE';
if robust==1
    robustStr='on';
    robustPlotStr='R';
elseif robust==0
    robustStr='off';
    robustPlotStr=[];
end

sigTestStatStr='Spearman';

minX=-4;
maxX=15;%xticks(-6:2:16);
minY=-5;
maxY=11;%yticks(-6:2:12);

%conts
b1=allBreaks(1);
b2=allBreaks(2);
b3=allBreaks(3);
b4=allBreaks(4);
b5=allBreaks(5);
b6=allBreaks(6);
breakpoints=[b2 b4];
if strcmp(regionStr,'vDLPFC')
    breakpoints=[b2];
end
regionBoundaries=[b1 b2;
                  b3+0.21 b4;
                  b5+0.21 b6];
jitterBoundary=1.5;

    %% Get RefAICc only
    if strcmp(regionStr,'dDLPFC')
        twoRegionElectrodes=find(x<=b4);
        x2region=x(twoRegionElectrodes);
        y2region=y(twoRegionElectrodes);
    elseif strcmp(regionStr,'vDLPFC')
        x2region=x;
        y2region=y;
    end

    %calculate each break, and residual of the resultant 2 lines flanking it
    best2BreakPoints=find2Breaks(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    best1BreakPoint=find1Break(x2region,y2region,breakpoints,regionBoundaries,jitterBoundary,robustStr);
    %fprintf('  2Breaks: %.1f and %.1f;    1Break: %.1f\n\n',best2BreakPoints(1),best2BreakPoints(end),best1BreakPoint)
    
    %% 1 segment model
    %fit model
    %twoRegionElectrodes=find(x<=b4);
    %x2region=x(twoRegionElectrodes);
    %y2region=y(twoRegionElectrodes);
    %x = x(:); % easier as column vectors
    %y = y(:);
    mdlA=fitlm(x2region,y2region,'RobustOpts',robustStr);
    switch statOfInterest
        case {'AIC'}
            %AIC
            SSR_1Seg = sum((mdlA.Residuals.Raw).^2);
            K_1Seg=mdlA.NumEstimatedCoefficients;
            N_1Seg=mdlA.NumObservations;
            AICc1=AIC(SSR_1Seg,N_1Seg,K_1Seg);
            %AICc1=999;
        case {'RMSE'}
            SSR_1Seg = sum((mdlA.Residuals.Raw).^2);
            K_1Seg=mdlA.NumEstimatedCoefficients;
            N_1Seg=mdlA.NumObservations;
            AICc1=AIC(SSR_1Seg,N_1Seg,K_1Seg);
            %AICc1=999;
            [arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);
    end
    
    
    %% 2 segment model
    %fit model
    xb1 =best1BreakPoint(1,1);
    if ~isnan(xb1)==0
        validBreak=0;
    else
        edges = [min(x2region),xb1,max(x2region)];
        bins = discretize(x2region,edges);
        coeffs=[];
        [idx1,~]=find(bins == 1);
        [idx2,~]=find(bins == 2);
        validBreak=~isempty(idx1)+~isempty(idx2)+sum(isnan(edges)); %check if empty and no-nan
    end
    
    switch validBreak
        case {2}
            %m1
            [idx,~]=find(bins == 1);
            x1=x2region(idx);
            y1=y2region(idx);
            mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
            %m2
            [idx,~]=find(bins == 2);
            x2=x2region(idx);
            y2=y2region(idx);
            mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
            %AIC
            switch statOfInterest
                case {'AIC'}
                    SSR_2Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2);
                    K_2Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients;
                    N_2Seg=mdlA.NumObservations+mdlB.NumObservations;
                    AICc2=AIC(SSR_2Seg,N_2Seg,K_2Seg);
                case {'RMSE'}
                    SSR_2Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2);
                    K_2Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients;
                    N_2Seg=mdlA.NumObservations+mdlB.NumObservations;
                    AICc2=AIC(SSR_2Seg,N_2Seg,K_2Seg);
                    [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
            end
        otherwise
            AICc2=NaN;
            arsq2whole=NaN;
    end

    %% 3 segment model
    %fit model
    x = x(:); % easier as column vectors
    y = y(:);
    xb1 =min(best2BreakPoints);
    xb2 =max(best2BreakPoints);
    
    if ~isnan(xb1)==0 || ~isnan(xb1)==0
        validBreak=0;
    else
        edges = [min(x),xb1,xb2,max(x)];
        bins = discretize(x,edges);
        coeffs=[];
        [idx1,~]=find(bins == 1);
        [idx2,~]=find(bins == 2);
        [idx3,~]=find(bins == 3);
        validBreak=~isempty(idx1)+~isempty(idx2)+~isempty(idx3)+sum(isnan(edges));
    end
    
    switch validBreak
        case {3}
            %m1
            x1=x(idx1);
            y1=y(idx1);
            mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
            %m2
            x2=x(idx2);
            y2=y(idx2);
            mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
            %m3
            x3=x(idx3);
            y3=y(idx3);
            mdlC=fitlm(x3,y3,'RobustOpts',robustStr);
            switch statOfInterest
                case {'AIC'}
                    %AIC
                    SSR_3Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2) + sum((mdlC.Residuals.Raw).^2);
                    K_3Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients+mdlC.NumEstimatedCoefficients;
                    N_3Seg=mdlA.NumObservations+mdlB.NumObservations+mdlC.NumObservations;
                    AICc3=AIC(SSR_3Seg,N_3Seg,K_3Seg);
                case {'RMSE'}
                    %AIC
                    SSR_3Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2) + sum((mdlC.Residuals.Raw).^2);
                    K_3Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients+mdlC.NumEstimatedCoefficients;
                    N_3Seg=mdlA.NumObservations+mdlB.NumObservations+mdlC.NumObservations;
                    AICc3=AIC(SSR_3Seg,N_3Seg,K_3Seg);
                    %[R2Adj_3Seg_whole,R2Adj_3Seg_perseg,RMSE_3Seg_whole,RMSE_3Seg_perseg]=getMdlStats(mdlA,mdlB,mdlC);
                    [arsq3whole,~,~,~]=getMdlStats(mdlA,mdlB,mdlC);
            end
        otherwise
            AICc3=NaN;
            arsq3whole=NaN;
    end
    arsqCont=[arsq1whole,arsq2whole,arsq3whole];
    %% Define reference AIC
    switch statOfInterest
        case {'AIC'}
            %refAICc=min([AICc1,AICc2,AICc3]);
        case {'RMSE'}
            
                shuffledDist1=shuffledAdjR2Cont(:,1);
                shuffledDist1 = shuffledDist1(~isnan(shuffledDist1));
                shuffledDist1_5prctile=prctile(shuffledDist1,95);
                %[H1,P1]=ttest(shuffledDist1,arsq1whole,'Tail','left');
                P1=1-invprctile(shuffledDist1,arsq1whole)/100;
                shuffledDist2=shuffledAdjR2Cont(:,2);
                shuffledDist2 = shuffledDist2(~isnan(shuffledDist2));
                shuffledDist2_5prctile=prctile(shuffledDist2,95);
                %[H2,P2]=ttest(shuffledDist2,arsq2whole,'Tail','left');
                P2=1-invprctile(shuffledDist2,arsq2whole)/100;
                shuffledDist3=shuffledAdjR2Cont(:,3);
                shuffledDist3 = shuffledDist3(~isnan(shuffledDist3));
                shuffledDist3_5prctile=prctile(shuffledDist3,95);
                %[H3,P3]=ttest(shuffledDist3,arsq3whole,'Tail','left');
                P3=1-invprctile(shuffledDist3,arsq3whole)/100;
                
                
                fprintf('TEST: p=%.2f (%.2f) | p=%.2f (%.2f) | p=%.2f (%.2f)\n',P1,shuffledDist1_5prctile,...
                    P2,shuffledDist2_5prctile,P3,shuffledDist3_5prctile)
                
                wantedregionStrsNew={'xLPFC','dDLPFC','vDLPFC'};
                %{
                if strcmp(regionStr,'dDLPFC')
                    [h, crit_p, ~, adj_p123]=fdr_bh([P1, P2, P3],.05,'pdep','no');
                    P1=adj_p123(1);
                    P2=adj_p123(2);
                    P3=adj_p123(3);
                elseif strcmp(regionStr,'vDLPFC')
                    [h, crit_p, ~, adj_p12]=fdr_bh([P1, P2],.05,'pdep','no');
                    P1=adj_p12(1);
                    P2=adj_p12(2);
                end
                %}
                if sum(~isnan([P1, P2, P3]))==3
                    %[h, crit_p, ~, adj_p123]=fdr_bh([P1, P2, P3],.05,'pdep','no');
                    %P1=adj_p123(1);
                    %P2=adj_p123(2);
                    %P3=adj_p123(3);
                    
                elseif sum(~isnan([P1, P2, P3]))==2
                    %[h, crit_p, ~, adj_p12]=fdr_bh([P1, P2],.05,'pdep','no');
                    %P1=adj_p12(1);
                    %P2=adj_p12(2);
                else
                    %
                end
                pvalues=[P1, P2, P3];
                
                
                if P1<.05
                    %add P1
                    AICcont=[AICcont,AICc1];
                end
                if P2<.05
                    %add P1
                    AICcont=[AICcont,AICc2];
                end
                if P3<.05
                    %add P1
                    AICcont=[AICcont,AICc3];
                end
            
            refAICc=min(AICcont); %only for sig models
    end
end









