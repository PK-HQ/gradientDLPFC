function [boundaryCellCont,fitsPCell]=fitDiscontinuousModelAdjP_2model(xCoords,yVal,...
anatBounds,yAxLabel,yLimits,yTickInterval,funcMeasureStr,...
surrData1seg,surrData2seg,surrData3seg,...
surrData1segSaveFilename,surrData2segSaveFilename,surrData3segSaveFilename,...
dvDLPFCStr,subplotIdx,...
adjGlobalGradientPvals,robustFitting,nRegions,loadFilteredData)

if contains(surrData1segSaveFilename,'J_sLPFC_Disselpct')
end
% This tests if 1-seg and 2-seg models have significant gradient, and
% identifies the better fitting model to the data along the anterior-posterior x-axis

%% DEFINE BASIC PARAMS

% CRITICAL
% Robust linear fitting
if robustFitting==1 
    robustStr='on';
elseif robustFitting==0
    robustStr='off';
end
% Define anatomical boundaries for functional graduient model fitting
breakpoints=[anatBounds(2) anatBounds(4)];
regionBoundaries=[anatBounds(1) anatBounds(2); %start and end of 8A
                  anatBounds(2)+0.01 anatBounds(4); %start and end of 9/46
                  anatBounds(4)+0.01 anatBounds(6)]; %start and end of 46
rangeBoundary=1.5; % constrain the functional boundary of fitted model to be +- 1.5 of anatomical boundary
adjGlobalGradientPvals(3)=NaN; % 2 regions only (areas 8A & 9/46), remove 3rd region (area 46)
boundaryCellCont=cell(1,3); % define cell to store region boundaries


% COSMETIC
lineColor=[0 0 0]; % line color for line plotting
dotSize=100; % datapoint diameter 
% Region colors and color transparency/alpha values
facealphaval=0.5;
if strcmp(dvDLPFCStr,'dDLPFC')
    regionShade1=[47 103 174]/255;
    regionShade2=[210 36 40]/255;
    regionShade3=[230 187 32]/255;
elseif strcmp(dvDLPFCStr,'vDLPFC')
    breakpoints=regionBoundaries(4);
    facealphaval=.5;
    regionShade1=[97 71 150]/255;
    regionShade2=[50 140 59]/255;
    regionShade3=[1 1 1]; %white/blank, no electrodes here (area 46)
end

%% CALCULATE ADJUSTED R2 OF EACH MODEL (1/2 SEGMENT MODEL)
%% DETERMINE IF THERE IS A SIGNIFICANT GLOBAL GRADIENT (SIGNIFICANT SHUFFLE TEST)
%% DETERMINE WHICH IS THE BETTER MODEL (SIGNIFICANT MODEL-COMPARISON P-VALUE FROM SURROGATE DATA TEST)

% Define electrodes for 8A and 9/46
if strcmp(dvDLPFCStr,'dDLPFC')
    electrodes2Region=find(xCoords<regionBoundaries(5)); %use all electrodes up to 9/46 (no 46)
    x2Region=xCoords(electrodes2Region);
    y2Region=yVal(electrodes2Region);
elseif strcmp(dvDLPFCStr,'vDLPFC')
    electrodes2Region=find(xCoords<9999999); %use all electrodes up to 9/46 (no 46)
    x2Region=xCoords(electrodes2Region);
    y2Region=yVal(electrodes2Region);
end

% Calculate the best functional boundary/break for the segmented models, by minimizing the RMSE
% of the overall fit. Break points can only be +-1.5 of the anatomical boundary position
% (i.e. constrained to a plausible position)
best1BreakPoint=find1Break(x2Region,y2Region,breakpoints,regionBoundaries,rangeBoundary,robustStr); % 2 segment fit

%% FIT 2-SEG MODEL, CALCULATE ADJ. R2, TEST FOR SIGNIFICANT GRADIENT
%fit model (only 8A and 946 datapoints used)
x2Region=xCoords(electrodes2Region);
y2Region=yVal(electrodes2Region);
mdlA=fitlm(x2Region,y2Region,'RobustOpts',robustStr);
[arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);

%% FIT 2-SEG MODEL, CALCULATE ADJ. R2, TEST FOR SIGNIFICANT GRADIENT
%fit model
bestBreakpoint=best1BreakPoint(1,1);

%If breakpoints are valid, split datapoints into 3 parts, 1 for each segment according to
%best breakpoints
if ~isnan(bestBreakpoint)==0
    validBreak=0;
else    
    edges = [min(x2Region),bestBreakpoint,max(x2Region)];
    bins = discretize(x2Region,edges);
    [idx1,~]=find(bins == 1);
    [idx2,~]=find(bins == 2);
    validBreak=~isempty(idx1)+~isempty(idx2)+sum(isnan(edges)); 
end

switch validBreak %check that break matrix is full (not-empty and no-nans)
    case {2}
        %fit segment 1
        [idx,~]=find(bins == 1);
        x1=x2Region(idx);
        y1=y2Region(idx);
        mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
        %fit segment 2
        [idx,~]=find(bins == 2);
        x2=x2Region(idx);
        y2=y2Region(idx);
        mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
        [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
    otherwise
        arsq2whole=NaN;
end

%% IDENTIFY BEST MODEL
[fitsCell,fitsPCell,dist1segRef_usableN,dist1segRef_range,dist2segRef_usableN,dist2segRef_range,...
    dist3segRef_usableN,dist3segRef_range]=whosthebest_v2(y2Region,adjGlobalGradientPvals,[arsq1whole,arsq2whole,NaN],...
    surrData1seg,surrData2seg,surrData3seg,...
    surrData1segSaveFilename,surrData2segSaveFilename,surrData3segSaveFilename,funcMeasureStr,nRegions,loadFilteredData);
if contains(surrData1segSaveFilename,'J_sLPFC_Disselpct') %no surrogate test done
    fitsPCell={'n/a', 'n/a', 'n/a'};
end
%assign adjusted p-values for 1/2/3 segment fits, for downstream
%plotting
P1=adjGlobalGradientPvals(1);
P2=adjGlobalGradientPvals(2);

%% PLOT 1/2/3 SEGMENT FITS, AND MARK SIGNIFICANCE.
%% SIGNIFICANT GRADIENT: BOLD ADJ R2 TEXT
%% BETTER MODEL: BOLD pMC VALUE   
% Note: if twoRegionElectrodes is used as an index, only 1-seg and 2-seg
% fits will be plotted and compared (i.e. only over areas 8A and 9/46).
% This is to maintain comparability across monkeys 1 and 2

%% PLOT 1-SEG MODEL
% fit lin model
mdl=fitlm(x2Region,y2Region,'RobustOpts',robustStr);
coeffs=mdl.Coefficients.Estimate;
% open subplot grid
subpltIdx=subplotIdx(1);
subplot(2,3,subpltIdx);hold on
% define x and y axis lims and their ticks
ylim([yLimits(1),yLimits(2)])
addSkippedTicks(yLimits(1):yTickInterval:yLimits(2),'y')
xlim([-4 14.01])
addSkippedTicks([-4:2:16],'x')
ylims=ylim;
% add figure alphabet label
%text(-6.9,ylims(2)+0.125*(ylims(2)-ylims(1)),subpltLabel(subplotIdx(1)),'FontSize',50,'Fontweight','bold');
% add transparent color shading for each region
area([min(xlim) regionBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade1);
area([regionBoundaries(4) anatBounds(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade2);
area([anatBounds(4) max(xlim)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade3);
% add anatomical boundaries as dashed black lines
line([regionBoundaries(4) regionBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
line([anatBounds(4) anatBounds(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
eqnY=coeffs(1)+coeffs(2)*[x2Region(1)-0.12 x2Region' x2Region(end)-0.1];
% plot datapoints
scatter(xCoords,yVal,dotSize,[1 1 1],'filled','MarkerEdgeColor','k');
% plot model fits
plot([x2Region(1)-0.12 x2Region' x2Region(end)-0.1],eqnY,'Color',lineColor,'LineWidth',3.5);
% save breakpoints for plotting breakpoint plot. 
%bestBreakpoint=[];
periBreakPoints=numel(x2Region);
if strcmp(fitsCell{1},'+')==1
    % Format: {significance of fit; break point; datapoints on each side of breakpoint}
    boundaryCellCont{1}={1;bestBreakpoint;periBreakPoints};
else
    boundaryCellCont{1}={0;bestBreakpoint;periBreakPoints};
end

% add adjusted r2 value either in bold (significant in shuffle test) or
% normal font (nonsig)
if P1<.05 & dist1segRef_usableN==1000 & fitsPCell{1}<.05 & isnan(bestBreakpoint)==0
    text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
elseif P1<.05 & dist1segRef_usableN==1000 & fitsPCell{1}>.05 & isnan(bestBreakpoint)==0
    text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','Normal','Interpreter','tex');
elseif contains(surrData1segSaveFilename,'J_sLPFC_Disselpct')
    text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','Normal','Interpreter','tex');
elseif P1<.05 & dist1segRef_usableN<1000 & isnan(bestBreakpoint)==0
    text(5.75,ylims(2)-0.295*(ylims(2)-ylims(1)),['N<1000, ' sprintf('%.0f',dist1segRef_usableN) fitsCell{1}],'Color',[0 0 0],'FontSize',6,'FontWeight','bold');
end

if P1<.001
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq1whole,2)) '^{***}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');    
elseif P1>.001 & P1<.01
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq1whole,2)) '^{**}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');    
elseif P1>.01 & P1<.05
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq1whole,2)) '^{*}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
elseif P1>.05
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq1whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13);
else
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2 = e' sprintf('%.02f',round(arsq1whole,2)) '^' round(P1,3) ' '],'Color',[0 0 0],'FontSize',13);
end


if subpltIdx==1 || subpltIdx==4
    ylabel(yAxLabel)
end
upFontSize(22,0.025);
addSkippedTicks([-4:2:16],'x')
hold off;

%% PLOT 2-SEG MODEL
bestBreakpoint =best1BreakPoint;
if isnan(bestBreakpoint)==0
    edges = [min(x2Region),bestBreakpoint,max(x2Region)];
    bins = discretize(x2Region,edges);
    coeffs=[];

    %fit first segment
    [idx,~]=find(bins == 1);
    x1=x2Region(idx);
    y1=y2Region(idx);
    mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
    coeffs=[coeffs;mdlA.Coefficients.Estimate(1);mdlA.Coefficients.Estimate(2)];

    %fit second segment
    [idx,~]=find(bins == 2);
    x2=x2Region(idx);
    y2=y2Region(idx);
    mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
    coeffs=[coeffs;mdlB.Coefficients.Estimate(1);mdlB.Coefficients.Estimate(2)];
    
    %define y-values of first and second segment fits
    y1=coeffs(1)+[x1(1)-0.12 x1' x1(end)-0.1]*coeffs(2);
    y2=coeffs(3)+[x2(1)-0.12 x2' x2(end)-0.1]*coeffs(4);
        
    if min([numel(x1) numel(x2)])<3 
        %don't plot if <3 datapoints on either segment
    else
        subplot(2,3,subplotIdx(2));hold on

        %draw region boundaries
        ylim([yLimits(1),yLimits(2)])
        addSkippedTicks(yLimits(1):yTickInterval:yLimits(2),'y')
        xlim([-4 14.01])

        addSkippedTicks([-4:2:16],'x')
        ylims=ylim;
       
        %add alphabet fig label
        %text(-6.9,ylims(2)+0.125*(ylims(2)-ylims(1)),subpltLabel(subplotIdx(2)),'FontSize',50,'Fontweight','bold');
        
        %shading and anat boundaries
        area([min(xlim) regionBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade1);
        area([regionBoundaries(4) anatBounds(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade2);
        area([anatBounds(4) max(xlim)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade3);
        line([regionBoundaries(4) regionBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
        line([anatBounds(4) anatBounds(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
        line([-5 15],0,'color',[0.75 0.75 0.75],'lineStyle','-','lineWidth',1)

        %plot datapoints and lines
        scatter(xCoords,yVal,dotSize,[1 1 1],'filled','MarkerEdgeColor','k')
        plot([x1(1)-0.12 x1' x1(end)-0.1],y1,'Color',lineColor,'LineWidth',3.5)
        plot([x2(1)-0.12 x2' x2(end)-0.1],y2,'Color',lineColor,'LineWidth',3.5);
        upFontSize(22,0.025);

        periBreakPoints=[numel(x1),numel(x2)];
        %save breakpoints
        if strcmp(fitsCell{2},'+')==1
            boundaryCellCont{2}={1;bestBreakpoint;periBreakPoints};
        else
            boundaryCellCont{2}={0;bestBreakpoint;periBreakPoints};
        end
        
        
        % add adj r2 bold or not depending on sig shuffle test
        if P1<.05 & dist2segRef_usableN==1000 & fitsPCell{2}<.05
            text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{2})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
        elseif P1<.05 & dist2segRef_usableN==1000 & fitsPCell{2}>.05
            text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{2})],'Color',[0 0 0],'FontSize',6,'FontWeight','normal','Interpreter','tex');
        elseif contains(surrData1segSaveFilename,'J_sLPFC_Disselpct')
            text(-3.8,ylims(2)-0.11*(ylims(2)-ylims(1)),['{\it p}_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','Normal','Interpreter','tex');
        elseif P1<.05 & dist2segRef_usableN<1000
            text(5.75,ylims(2)-0.295*(ylims(2)-ylims(1)),['N<1000,' sprintf('%.0f',dist2segRef_usableN) fitsCell{2}],'Color',[0 0 0],'FontSize',6,'FontWeight','bold');
        end

        if P2<.001
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq2whole,2)) '^{***}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
        elseif P2>.001 & P2<.01
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq2whole,2)) '^{**}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
        elseif P2>.01 & P2<.05
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq2whole,2)) '^{*}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
        elseif P2>.05           
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2=' sprintf('%.02f',round(arsq2whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13);
        else
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. R^2 = e' sprintf('%.02f',round(arsq2whole,2))],'Color',[0 0 0],'FontSize',13);
        end
    end
end
upFontSize(22,0.025);
hold off

if min(yLimits)>min(y2Region) || max(yLimits)<max(y2Region)
    fprintf('WRONG Y-LIMS: IT SHOULD BE [%.2f %.2f] INSTEAD OF [%.2f %.2f]\n',min(y2Region),max(y2Region),min(yLimits),max(yLimits)) 
    1;
end
end