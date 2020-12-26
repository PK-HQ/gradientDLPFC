function [breakpointCellCont,fitsPCell]=fitDiscontinuousModelAdjP_3model(x,y,...
anatomicalBoundaries,yLab,statLimits,statSpacing,statName,...
surrogateAdjR2Cont1seg,surrogateAdjR2Cont2seg,surrogateAdjR2Cont3seg,...
surrogate1segFname,surrogate2segFname,surrogate3segFname,...
regionStr,subplotIdx,...
adjP,robust,nRegionStr,loadPrefilteredData)
regionColor=[0 0 0];
figure
%% Define some params
breakpointCellCont=cell(1,3); % define breakpoint storage cell
dotSize=100; %datapoint diameter 
if robust==1
    robustStr='on';
elseif robust==0
    robustStr='off';
end

% define anatomical breakpoints
breakpoints=[anatomicalBoundaries(2) anatomicalBoundaries(4)];
regionBoundaries=[anatomicalBoundaries(1) anatomicalBoundaries(2);
                  anatomicalBoundaries(2)+0.01 anatomicalBoundaries(4);
                  anatomicalBoundaries(4)+0.01 anatomicalBoundaries(6)];
jitterBoundary=1.5; %boundary +- 1.5 is where the fitted model's functional breakpoint can exist

% def region shading colors and transparency
facealphaval=0.5;
if strcmp(regionStr,'sLPFC')
    regionShade1=[47 103 174]/255;
    regionShade2=[210 36 40]/255;
    regionShade3=[230 187 32]/255;
elseif strcmp(regionStr,'iLPFC')
    breakpoints=regionBoundaries(4);
    facealphaval=.5;
    regionShade1=[97 71 150]/255;
    regionShade2=[50 140 59]/255;
    regionShade3=[1 1 1]; %white/blank, no electrodes here
end

%% CALCULATE ADJUSTED R2 OF EACH MODEL (1/2/3 SEGMENTS), AND DETERMINE THE BEST FITTING MODEL 
%% USING THE SHUFFLE TEST + SURROGATE DATA TEST

% calculate the best functional boundary/break for the segmented models, by minimizing the RMSE of the overall fit
% break points can only be +-1.5 of the anatomical boundary position
% (i.e. constrained to a plausible position)

if strcmp(regionStr,'sLPFC')
    allRegionElectrodes=find(x<99999);
    x2region=x(allRegionElectrodes);
    y2region=y(allRegionElectrodes);
    %adjP(3)=NaN;
elseif strcmp(regionStr,'iLPFC')
    allRegionElectrodes=find(x<99999); %all
    x2region=x(allRegionElectrodes);
    y2region=y(allRegionElectrodes);
end

[best2BreakPoints,tiedBreakpoints2Breaks]=find2Breaks(x,y,breakpoints,regionBoundaries,jitterBoundary,robustStr); % 3 segment fit (unused)
[best1BreakPoint,tiedBreakpoints1Breaks]=find1Break(x2region,y2region,breakpoints,regionBoundaries,jitterBoundary,robustStr); % 2 segment fit

%% Fit 1 segment model and calculate Adj. R2
%fit model (only 8A and 946 datapoints used)
x2region=x(allRegionElectrodes);
y2region=y(allRegionElectrodes);
mdlA=fitlm(x2region,y2region,'RobustOpts',robustStr);
[arsq1whole,~,~,~]=getMdlStats(mdlA,[],[]);

%% Fit 2 segment model and calculate Adj. R2
%fit model
%x2region=x(twoRegionElectrodes);
%y2region=y(twoRegionElectrodes);
xb1 =best1BreakPoint(1,1);

%If breakpoints are valid, split datapoints into 3 parts, 1 for each segment according to
%best breakpoints
if ~isnan(xb1)==0
    validBreak=0;
else    
    edges = [min(x2region),xb1,max(x2region)];
    bins = discretize(x2region,edges);
    [idx1,~]=find(bins == 1);
    [idx2,~]=find(bins == 2);
    validBreak=~isempty(idx1)+~isempty(idx2)+sum(isnan(edges)); 
end

switch validBreak %check that break matrix is full (not-empty and no-nans)
    case {2}
        %fit segment 1
        [idx,~]=find(bins == 1);
        x1=x2region(idx);
        y1=y2region(idx);
        mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
        %fit segment 2
        [idx,~]=find(bins == 2);
        x2=x2region(idx);
        y2=y2region(idx);
        mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
        [arsq2whole,~,~,~]=getMdlStats(mdlA,mdlB,[]);
    otherwise
        arsq2whole=NaN;
end

%% Fit 3 segment model and calculate Adj. R2
%fit model
%x2region = x(:); % easier as column vectors
%y2region = y(:);
xb1 =min(best2BreakPoints);
xb2=max(best2BreakPoints);

%If breakpoints are valid, split datapoints into 3 parts, 1 for each segment according to
%best breakpoints
if ~isnan(xb1)==0 || ~isnan(xb1)==0
    validBreak=0;
else
    edges = [min(x),xb1,xb2,max(x)];
    bins = discretize(x,edges);
    [idx1,~]=find(bins == 1);
    [idx2,~]=find(bins == 2);
    [idx3,~]=find(bins == 3);
    validBreak=~isempty(idx1)+~isempty(idx2)+~isempty(idx3)+sum(isnan(edges));
end

switch validBreak %check that break matrix is full (not-empty and no-nans)
    case {3}
        %fit segment 1
        x1=x(idx1);
        y1=y(idx1);
        mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
        %fit segment 2
        x2=x(idx2);
        y2=y(idx2);
        mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
        %fit segment 3
        x3=x(idx3);
        y3=y(idx3);
        mdlC=fitlm(x3,y3,'RobustOpts',robustStr);
        [arsq3whole,~,~,~]=getMdlStats(mdlA,mdlB,mdlC);
    otherwise
        arsq3whole=NaN;
end

%% Find the best model with surrogate test
%{
[fitsCell,dist1segRef_usableN,dist1segRef_range,dist2segRef_usableN,dist2segRef_range,...
    dist3segRef_usableN,dist3segRef_range]=whosthebest(y2region,adjP,[arsq1whole,arsq2whole,arsq3whole],...
    surrogateAdjR2Cont1seg,surrogateAdjR2Cont2seg,surrogateAdjR2Cont3seg,statName);
%}
[fitsCell,fitsPCell,dist1segRef_usableN,dist1segRef_range,dist2segRef_usableN,dist2segRef_range,...
    dist3segRef_usableN,dist3segRef_range]=whosthebest_v2(y2region,adjP,[arsq1whole,arsq2whole,arsq3whole],...
    surrogateAdjR2Cont1seg,surrogateAdjR2Cont2seg,surrogateAdjR2Cont3seg,...
    surrogate1segFname,surrogate2segFname,surrogate3segFname,statName,nRegionStr,loadPrefilteredData);

%{
[fitsCell,dist1segRef_usableN,dist1segRef_range,dist2segRef_usableN,dist2segRef_range,...
    dist3segRef_usableN,dist3segRef_range]=whosthebest(y2region,adjP,[arsq1whole,arsq2whole,arsq3whole],...
    surrogateAdjR2Cont1seg,surrogateAdjR2Cont2seg,surrogateAdjR2Cont3seg,...
    surrogate1segFname,surrogate2segFname,surrogate3segFname,statName);
    %}
%assign adjusted p-values for 1/2/3 segment fits, for downstream
%plotting
P1=adjP(1);
P2=adjP(2);
P3=adjP(3);




%% PLOT 1/2/3 SEGMENT FITS, AND MARK SIGNIFICANCE.
%% SIG. SHUFFLE: BOLD ADJ R2 TEXT vs NORMAL FONT
%% SIG. SURROGATE (BEST MODEL): THICKER BLACK LINE vs GREY LINE    
% Note: if twoRegionElectrodes is used as an index, only linear and 2 segment
% fits will be plotted, and only over areas 8A and 9/46. This is to
% maintain comparability across monkeys 1 and 2

%% 1 segment model
% define x, y datapoints used for plotting in 8A and 9/46
%twoRegionElectrodes=find(x<=anatomicalBoundaries(4));
x2region=x(allRegionElectrodes);
y2region=y(allRegionElectrodes);
%x = x(:); % easier as column vectors
%y = y(:);

% fit lin model
mdl=fitlm(x2region,y2region,'RobustOpts',robustStr);
coeffs=mdl.Coefficients.Estimate;
% open subplot grid
subpltIdx=subplotIdx(1);
subplot(2,3,subpltIdx);hold on
% define x and y axis lims and their ticks
ylim([statLimits(1),statLimits(2)])
addSkippedTicks(statLimits(1):statSpacing:statLimits(2),'y')
xlim([-4 14.01])
addSkippedTicks([-4:2:16],'x')
ylims=ylim;
% add figure alphabet label
text(-6.9,ylims(2)+0.125*(ylims(2)-ylims(1)),subpltLabel(subplotIdx(1)),'FontSize',50,'Fontweight','bold');
% add transparent color shading for each region
area([min(xlim) regionBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade1);
area([regionBoundaries(4) anatomicalBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade2);
area([anatomicalBoundaries(4) max(xlim)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade3);
% add anatomical boundaries as dashed black lines
line([regionBoundaries(4) regionBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
line([anatomicalBoundaries(4) anatomicalBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
eqnY=coeffs(1)+coeffs(2)*[x2region(1)-0.12 x2region' x2region(end)-0.1];
textColor=[0 0 0];
dataPointColor=[1 1 1];
lineColor=[0 0 0];
% mark lines with black or grey depending on whether its the best model or
% not
if P1<.05
    lineColor=[0 0 0];
elseif P1>.05
    lineColor=[0 0 0];
end
%{
if P1<.05
    if strcmp(fitsCell{2},'+')==0 & P2<.05 & strcmp(fitsCell{3},'+')==0  & P3<.05 || P2>.05 & P3>.05
        lineColor=[0 0 0];
    elseif strcmp(fitsCell{1},'+')==1
        lineColor=[0 0 0];
    elseif strcmp(fitsCell{1},'+')==0
        lineColor=[0 0 0];
    end
end
%}

% plot datapoints
scatter(x,y,dotSize,[1 1 1],'filled','MarkerEdgeColor','k');
% plot model fits
plot([x2region(1)-0.12 x2region' x2region(end)-0.1],eqnY,'Color',lineColor,'LineWidth',3.5);
% save breakpoints for plotting breakpoint plot. 
% Format: {significance of fit; break point; datapoints on each side of breakpoint}
breakpos=[];
periBreakPoints=numel(x2region);

if P1<.05%strcmp(fitsCell{1},'+')==1
    breakpointCellCont{1}={1;breakpos;periBreakPoints;tiedBreakpoints1Breaks};
else
    breakpointCellCont{1}={0;breakpos;periBreakPoints};
end

% add adjusted r2 value either in bold (significant in shuffle test) or
% normal font (nonsig)
p1v2and3=fitsPCell{1};
p1v2=p1v2and3(1);
p1v3=p1v2and3(2);
if P1<.05 & dist1segRef_usableN==1000 & p1v2<0.05 & p1v3<0.05
    text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 1-v-2}=' sprintf('%.2f',p1v2)],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
    text(-3.6,ylims(2)-0.255*(ylims(2)-ylims(1)),['p_{MC 1-v-3}=' sprintf('%.2f',p1v3)],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
elseif P1<.05 & dist1segRef_usableN==1000 & p1v2>0.05 || P1<.05 & dist1segRef_usableN==1000 & p1v3>0.05
    text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 1-v-2}=' sprintf('%.2f',p1v2)],'Color',[0 0 0],'FontSize',6,'FontWeight','normal','Interpreter','tex');
    text(-3.6,ylims(2)-0.255*(ylims(2)-ylims(1)),['p_{MC 1-v-3}=' sprintf('%.2f',p1v3)],'Color',[0 0 0],'FontSize',6,'FontWeight','normal','Interpreter','tex');
elseif P1<.05 & dist1segRef_usableN<1000 & isnan(xb1)==0
    text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),'N<1000','Color',[0 0 0],'FontSize',6);
    text(-3.6,ylims(2)-0.255*(ylims(2)-ylims(1)),['p_{MC 1-v-2}=' sprintf('%.2f',p1v2)],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
    text(-3.6,ylims(2)-0.305*(ylims(2)-ylims(1)),['p_{MC 1-v-3}=' sprintf('%.2f',p1v3)],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
end
if P1<.05 & dist1segRef_usableN<1000
    %text(5.75,ylims(2)-0.095*(ylims(2)-ylims(1)),[sprintf('%.0f',dist1segRef_usableN) fitsCell{1}],'Color',[0 0 0],'FontSize',6,'FontWeight','bold');
end
if P1<.001
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq1whole,2)) '^{***}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');    
    
    %add text for size of distribution for this model, and min max R2
    %of distribution
    %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist1segRef_usableN,2)) ' (' sprintf('%.02f',round(dist1segRef_range(1),2)) ',' sprintf('%.02f',round(dist1segRef_range(2),2)) ')'],'FontSize',20,'Color',color1seg); 
elseif P1>.001 & P1<.01
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq1whole,2)) '^{**}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');    
    
    %add text for size of distribution for this model, and min max R2
    %of distribution
    %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist1segRef_usableN,2)) ' (' sprintf('%.02f',round(dist1segRef_range(1),2)) ',' sprintf('%.02f',round(dist1segRef_range(2),2)) ')'],'FontSize',20,'Color',color1seg); 
elseif P1>.01 & P1<.05
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq1whole,2)) '^{*}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
    
    %add text for size of distribution for this model, and min max R2
    %of distribution
    %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist1segRef_usableN,2)) ' (' sprintf('%.02f',round(dist1segRef_range(1),2)) ',' sprintf('%.02f',round(dist1segRef_range(2),2)) ')'],'FontSize',20,'Color',color1seg); 
elseif P1>.05
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq1whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13);
else
    text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2 = e' sprintf('%.02f',round(arsq1whole,2)) '^' round(P1,3) ' '],'Color',[0 0 0],'FontSize',13);
end


if subpltIdx==1 || subpltIdx==4
    ylabel(yLab)
end
upFontSize(22,0.025);
addSkippedTicks([-4:2:16],'x')
hold off;

%% 2 segment model
%fit model
%twoRegionElectrodes=find(x<=anatomicalBoundaries(4));
%x2region=x(twoRegionElectrodes);
%y2region=y(twoRegionElectrodes);
x2region=x(allRegionElectrodes);
y2region=y(allRegionElectrodes);


xb1 =best1BreakPoint;
if isnan(xb1)==0
    edges = [min(x2region),xb1,max(x2region)];
    bins = discretize(x2region,edges);
    %M = [x2region.*(bins == 1),bins == 1,x2region.*(bins == 2),bins == 2];

    coeffs=[];

    %fit first segment
    [idx,~]=find(bins == 1);
    x1=x2region(idx);
    y1=y2region(idx);
    mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
    coeffs=[coeffs;mdlA.Coefficients.Estimate(1);mdlA.Coefficients.Estimate(2)];

    %fit second segment
    [idx,~]=find(bins == 2);
    x2=x2region(idx);
    y2=y2region(idx);
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
        ylim([statLimits(1),statLimits(2)])
        addSkippedTicks(statLimits(1):statSpacing:statLimits(2),'y')
        xlim([-4 14.01])

        addSkippedTicks([-4:2:16],'x')
        ylims=ylim;
        xlims=xlim;
        
        %add alphabet fig label
        text(-6.9,ylims(2)+0.125*(ylims(2)-ylims(1)),subpltLabel(subplotIdx(2)),'FontSize',50,'Fontweight','bold');
        
        %shading and anat boundaries
        area([min(xlim) regionBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade1);
        area([regionBoundaries(4) anatomicalBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade2);
        area([anatomicalBoundaries(4) max(xlim)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade3);
        line([regionBoundaries(4) regionBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
        line([anatomicalBoundaries(4) anatomicalBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
        line([-5 15],0,'color',[0.75 0.75 0.75],'lineStyle','-','lineWidth',1)

        lineColor=[0 0 0];
        if P2<.05
            lineColor=[0 0 0];
        elseif P2>.05
            lineColor=[0 0 0];
        end
        %{
        if P2<.05
            if strcmp(fitsCell{2},'+')==1
                lineColor=[0 0 0];
            elseif strcmp(fitsCell{2},'+')==0
                lineColor=[0 0 0];
            end
        end
        %}

        %plot datapoints and lines
        scatter(x,y,dotSize,[1 1 1],'filled','MarkerEdgeColor','k')
        plot([x1(1)-0.12 x1' x1(end)-0.1],y1,'Color',lineColor,'LineWidth',3.5)
        plot([x2(1)-0.12 x2' x2(end)-0.1],y2,'Color',lineColor,'LineWidth',3.5);
        upFontSize(22,0.025);

        breakpos=xb1;
        periBreakPoints=[numel(x1),numel(x2)];
        
        %save breakpoints
        if P2<.05%strcmp(fitsCell{2},'+')==1
            breakpointCellCont{2}={1;breakpos;periBreakPoints;tiedBreakpoints1Breaks};
        else
            breakpointCellCont{2}={0;breakpos;periBreakPoints;tiedBreakpoints1Breaks};
        end
        
        if P1<.05 & dist2segRef_usableN==1000 & fitsPCell{2}<0.05
            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 2-v-1}=' sprintf('%.2f',fitsPCell{2})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
        elseif P1<.05 & dist2segRef_usableN==1000 & fitsPCell{2}>0.05
            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 2-v-1}=' sprintf('%.2f',fitsPCell{2})],'Color',[0 0 0],'FontSize',6,'FontWeight','normal','Interpreter','tex');
            %text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
        elseif P1<.05 & dist2segRef_usableN<1000
            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),'N<1000','Color',[0 0 0],'FontSize',6);
            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 2-v-1}=' sprintf('%.2f',fitsPCell{2})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
        end
        if P1<.05 & dist2segRef_usableN<1000
            %text(5.75,ylims(2)-0.095*(ylims(2)-ylims(1)),[sprintf('%.0f',dist2segRef_usableN) fitsCell{2}],'Color',[0 0 0],'FontSize',6,'FontWeight','bold');
        end
        % add adj r2 bold or not depending on sig shuffle test
        if P2<.001
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq2whole,2)) '^{***}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
            
            %add text for size of distribution for this model, and min max R2
            %of distribution
            %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist2segRef_usableN,2)) ' (' sprintf('%.02f',round(dist2segRef_range(1),2)) ',' sprintf('%.02f',round(dist2segRef_range(2),2)) ')'],'FontSize',20,'Color',color2seg); 
        elseif P2>.001 & P2<.01
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq2whole,2)) '^{**}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
            
            %add text for size of distribution for this model, and min max R2
            %of distribution
            %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist2segRef_usableN,2)) ' (' sprintf('%.02f',round(dist2segRef_range(1),2)) ',' sprintf('%.02f',round(dist2segRef_range(2),2)) ')'],'FontSize',20,'Color',color2seg); 
        elseif P2>.01 & P2<.05
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq2whole,2)) '^{*}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
            
            %add text for size of distribution for this model, and min max R2
            %of distribution
            %text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist2segRef_usableN,2)) ' (' sprintf('%.02f',round(dist2segRef_range(1),2)) ',' sprintf('%.02f',round(dist2segRef_range(2),2)) ')'],'FontSize',20,'Color',color2seg); 
        elseif P2>.05           
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq2whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13);
        else
            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2 = e' sprintf('%.02f',round(arsq2whole,2))],'Color',[0 0 0],'FontSize',13);
        end
    end
end
upFontSize(22,0.025);
hold off
1;
%% 3 segment model (not used when only the posterior 2 regions of 8A and 9/46 are plotted) 
switch strcmp(regionStr,'dL')
    case {1}
        %skip this
    otherwise 
        %fit model
        x = x(:); % easier as column vectors
        y = y(:);
        xb1 =min(best2BreakPoints);
        xb2 =max(best2BreakPoints);
        if isnan(xb1)==0 & isnan(xb2)==0
            edges = [min(x),xb1,xb2,max(x)];
            bins = discretize(x,edges);
            M = [x.*(bins == 1),x.*(bins == 2),x.*(bins == 3)];

            coeffs=[];
            [idx1,~]=find(bins == 1);
            [idx2,~]=find(bins == 2);
            [idx3,~]=find(bins == 3);


            validBreak=~isempty(idx1)+~isempty(idx2)+~isempty(idx3);
            switch validBreak
                case {3}
                    %m1
                    x1=x(idx1);
                    y1=y(idx1);
                    mdlA=fitlm(x1,y1,'RobustOpts',robustStr);
                    coeffs=[coeffs;mdlA.Coefficients.Estimate(1);mdlA.Coefficients.Estimate(2)];
                    %[R1,P1] = corrcoef(x1,y1);

                    %m2
                    x2=x(idx2);
                    y2=y(idx2);
                    mdlB=fitlm(x2,y2,'RobustOpts',robustStr);
                    coeffs=[coeffs;mdlB.Coefficients.Estimate(1);mdlB.Coefficients.Estimate(2)];
                    %[R2,P2] = corrcoef(x2,y2);

                    %m3
                    x3=x(idx3);
                    y3=y(idx3);
                    mdlC=fitlm(x3,y3,'RobustOpts',robustStr);
                    coeffs=[coeffs;mdlC.Coefficients.Estimate(1);mdlC.Coefficients.Estimate(2)];
                    %[R3,P3] = corrcoef(x3,y3);

                    %SSR_3Seg = sum((mdlA.Residuals.Raw).^2) + sum((mdlB.Residuals.Raw).^2) + sum((mdlC.Residuals.Raw).^2);
                    %K_3Seg=mdlA.NumEstimatedCoefficients+mdlB.NumEstimatedCoefficients+mdlC.NumEstimatedCoefficients;
                    %N_3Seg=mdlA.NumObservations+mdlB.NumObservations+mdlC.NumObservations;


                    %AICc3=AIC(SSR_3Seg,N_3Seg,K_3Seg);


                    y1=coeffs(1)+[x1(1)-0.12 x1' x1(end)-0.1]*coeffs(2);

                    y2=coeffs(3)+[x2(1)-0.12 x2' x2(end)-0.1]*coeffs(4);

                    y3=coeffs(5)+[x3(1)-0.12 x3' x3(end)-0.1]*coeffs(6);
                    if min([numel(x1) numel(x2) numel(x3)])<3
                        %dont plot
                    elseif strcmp(regionStr,'iLPFC')==1
                        fprintf('istrue\n')
                    else

                        %figure;hold on
                        h1=subplot(2,3,subplotIdx(3));hold on

                        %draw region boundaries
                        %ylim([ylimits(1),ylimits(2)])
                        %addSkippedTicks(ylimits(1):tickInt:ylimits(2),'y')
                        ylim([statLimits(1),statLimits(2)])
                        addSkippedTicks(statLimits(1):statSpacing:statLimits(2),'y')
                        xlim([-4 14.01])
                        %xticks([minX:2:maxX])
                        %names = {'-8'; '-4'; '0'; '4'; '8'; '12'; '16'};
                        %set(gca,'xtick',[-8:4:16],'xticklabel',names)
                        addSkippedTicks([-4:2:16],'x')
                        xlims=xlim;
                        ylims=ylim;
                        %add alphabet fig label
                        text(-6.9,ylims(2)+0.125*(ylims(2)-ylims(1)),subpltLabel(subplotIdx(3)),'FontSize',50,'Fontweight','bold');
                        area([min(xlim) regionBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade1);
                        area([regionBoundaries(4) anatomicalBoundaries(4)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade2);
                        area([anatomicalBoundaries(4) max(xlim)],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',facealphaval, 'EdgeAlpha',0, 'FaceColor',regionShade3);
                        line([regionBoundaries(4) regionBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
                        line([anatomicalBoundaries(4) anatomicalBoundaries(4)],ylim,'color',[1 1 1],'lineStyle','--','lineWidth',2.5)
                        line([-5 15],0,'color',[0.75 0.75 0.75],'lineStyle','-','lineWidth',1)

                        textColor='black';
                        dataPointColor=[1 1 1];
                        lineColor=[0 0 0];
                        if P3<.05
                            lineColor=[0 0 0];
                        elseif P3>.05
                            lineColor=[0 0 0];
                        end
                        %{
                        if P3>.05
                            regionColor3=[0.75 0.75 0.75];
                        else
                            regionColor3=regionColor;
                            if strcmp(fitsCell{3},'+')==1
                                textColor=[0 0 0];
                                dataPointColor=[1 1 1];
                                lineColor=[0 0 0];
                            elseif strcmp(fitsCell{3},'+')==0
                                dataPointColor=[1 1 1];
                                lineColor=[0 0 0];
                            end
                        end
                        %}

                        scatter(x,y,dotSize,[1 1 1],'filled','MarkerEdgeColor','k');
                        plot([x1(1)-0.12 x1' x1(end)-0.1],y1,'Color',lineColor,'LineWidth',3.5);
                        plot([x2(1)-0.12 x2' x2(end)-0.1],y2,'Color',lineColor,'LineWidth',3.5);
                        plot([x3(1)-0.12 x3' x3(end)-0.1],y3,'Color',lineColor,'LineWidth',3.5);
                        upFontSize(22,0.025);

                        breakpos=[xb1,xb2];
                        periBreakPoints=[numel(x1),numel(x2),numel(x3)];
                        %{
                        if AICc3-refAICc<10
                            breakpointCellCont{3}={1;breakpos;periBreakaPoints};
                        else
                            breakpointCellCont{3}={0;breakpos;periBreakPoints};
                        end
                        %}

                        %[H3,P3]=ttest(abs(shuffledDist),abs(arsq3whole),'Tail','left');
                        %[H3,P3]=ttest(arsq3whole,shuffledAdjR2Cont(:,3),'Tail','right');
                        %%text(-3.9,ylims(2)-0.25*(ylims(2)-ylims(1)),[sprintf('%.02f',round(dist3segRef_usableN,2)) ' (' sprintf('%.02f',round(dist3segRef_range(1),2)) ',' sprintf('%.02f',round(dist3segRef_range(2),2)) ')'],'FontSize',20,'Color',[0 0 0]); 
                        %save breakpoints
                        if P3<.05%strcmp(fitsCell{3},'+')==1
                            breakpointCellCont{3}={1;breakpos;periBreakPoints;tiedBreakpoints2Breaks};
                        else
                            breakpointCellCont{3}={0;breakpos;periBreakPoints;tiedBreakpoints2Breaks};
                        end
                        if P1<.05 & dist3segRef_usableN==1000 & fitsPCell{3}<0.05
                            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 3-v-1}=' sprintf('%.2f',fitsPCell{3})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
                            %text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC}=' sprintf('%.2f',fitsPCell{1})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
                        elseif P1<.05 & dist2segRef_usableN==1000 & fitsPCell{3}>0.05
                            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 3-v-1}=' sprintf('%.2f',fitsPCell{3})],'Color',[0 0 0],'FontSize',6,'FontWeight','normal','Interpreter','tex');
                        elseif P1<.05 & dist3segRef_usableN<1000
                            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),'N<1000','Color',[0 0 0],'FontSize',6);
                            text(-3.6,ylims(2)-0.105*(ylims(2)-ylims(1)),['p_{MC 3-v-1}=' sprintf('%.2f',fitsPCell{3})],'Color',[0 0 0],'FontSize',6,'FontWeight','bold','Interpreter','tex');
                        end
                        if P1<.05 & dist3segRef_usableN<1000
                            %text(5.75,ylims(2)-0.095*(ylims(2)-ylims(1)),[sprintf('%.0f',dist3segRef_usableN) fitsCell{3}],'Color',[0 0 0],'FontSize',6,'FontWeight','bold');
                        end
                        if P3<.001
                            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq3whole,2)) '^{***}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            %text(-1.65,ylims(2)-0.2*(ylims(2)-ylims(1)),['{\Delta}BIC  = ' sprintf('%.0f',round(AICc3-refAICc,0))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            sig=1;
                            
                        elseif P3>.001 & P3<.01
                            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq3whole,2)) '^{**}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            %text(-1.65,ylims(2)-0.2*(ylims(2)-ylims(1)),['{\Delta}BIC  = ' sprintf('%.0f',round(AICc3-refAICc,0))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            sig=1;
                            breakpos=[xb1,xb2];
                            periBreakPoints=[numel(x1),numel(x2),numel(x3)];
                            
                        elseif P3>.01 & P3<.05
                            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq3whole,2)) '^{*}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            %text(-1.65,ylims(2)-0.2*(ylims(2)-ylims(1)),['{\Delta}BIC  = ' sprintf('%.0f',round(AICc3-refAICc,0))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            sig=1;
                            breakpos=[xb1,xb2];
                            periBreakPoints=[numel(x1),numel(x2),numel(x3)];
                            
                        elseif P3>.05
                            %text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq3whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(arsq3whole,2)) '^{n.s.}'],'Color',[0 0 0],'FontSize',13);
                            sig=0;
                            breakpos=[xb1,xb2];
                            periBreakPoints=[numel(x1),numel(x2),numel(x3)];
                            
                        else
                            text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2 = error' sprintf('%.02f',round(arsq3whole,2))],'Color',[0 0 0],'FontSize',13);
                        end
                    end



                    %breakpointCellCont{3}={sig;breakpos;periBreakPoints};
                    %text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(RMSE3Seg,2))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');

                    %text(-1.5,ylims(2)-0.16*(ylims(2)-ylims(1)),['{\Delta}_{\iti} = ' sprintf('%.02f',round(AICc3-refAICc,0))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');

                    %%text(xlims(1)+0.75*(xlims(2)-xlims(1)),ylims(1)+0.75*(ylims(2)-ylims(1)),['R = ' sprintf('%.02f',round(R1,2)) ' (p = ' sprintf('%.02f',round(P1(2),2)) ')'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                    %text(xlims(1)+0.20*(xlims(2)-xlims(1)),ylims(1)+0.75*(ylims(2)-ylims(1)),['R = ' sprintf('%.02f',round(R2,2)) ' (p = ' sprintf('%.02f',round(P2,2)) ')'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');
                    %text(xlims(1)+0.25*(xlims(2)-xlims(1)),ylims(1)+0.75*(ylims(2)-ylims(1)),['R = ' sprintf('%.02f',round(R3,2)) ' (p = ' sprintf('%.02f',round(P3,2)) ')'],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');

                    %text(-3.8,ylims(2)+0.095*(ylims(2)-ylims(1)),['Adj. r^2=' sprintf('%.02f',round(slm.stats.R2Adj,2))],'Color',[0 0 0],'FontSize',13,'FontWeight','bold');

                    %line([xb1 xb1],ylim,'color','r','lineStyle','-','lineWidth',.5)
                    %line([xb2 xb2],ylim,'color','r','lineStyle','-','lineWidth',.5)
                    %ylabel(yLab);xlabel(xLab);
                    xlim([-4 14.01])
                    %xticks([minX:2:maxX])
                    %names = {'-8'; '-4'; '0'; '4'; '8'; '12'; '16'};
                    %set(gca,'xtick',[-8:4:16],'xticklabel',names)
                    addSkippedTicks([-4:2:16],'x')

                    %title(wantedregion2,'fontweight','normal')
                    upFontSize(22,0.025);
                    hold off
                    upFontSize(22,0.025);

                    if subplotIdx(3)==6
                        %delete(h1)
                    end
                otherwise
                    %
            end
        end
end
%}

if min(statLimits)>min(y2region) || max(statLimits)<max(y2region)
    fprintf('ALERTALERTALERTALERTALERTALERTALERTALERT [%.2f %.2f] instead of [%.2f %.2f]\n',min(y2region),max(y2region),min(statLimits),max(statLimits)) 
    1;
end
end









