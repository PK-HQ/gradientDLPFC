function [a8AD,a946D,a8AV,a946V]=getFunctionalBoundsIndiv(Monky)
jitterBoundary=1.5;
statType='r';
order=fliplr([2 3 17 9 16 18 4 5 8 10 11 20 21]);%[17 9 10 11 18 16 20 21 14 15 12];
%surr_J__rbreakpointMeasuresCont_p005.mat
if Monky==1
    MonkyID='M2_';
    RegionsUsed=1; %only sDLPFC, as iDLPFC is only 8Av
    
    sDLPFCstatsUsed=order;
    iDLPFCstatsUsed=order;

    %% Load data and set params
    % Load single monkey data
    %SurrV2_' nRegionStr 'regions_' MonkyID '_rbreakpointMeasuresCont
    load(['~/Desktop/HPC_PK/data/em/3regions_' MonkyID 'breakpoints.mat'])
    %load(['~/Desktop/HPC_PK/data/em/SurrV2_3regions_' MonkyID '_rbreakpointMeasuresCont_p005.mat'])
    %load(['~/Desktop/HPC_PK/data/em/3Bregion_surr_' MonkyID '_' statType 'breakpointMeasuresCont_p005.mat'])
    %
    oneBreakCont=[];
    twoBreakCont=[];
    a8AD=999;
    a946D=10;
    a8AV=99;
    a946V=99;
    statNames={'All','Tarselpct','Disselpct','D1selpct','D2selpct','CS','LMS','NMSpct',...
        'Trf','D1rf','D2rf','Mixedstr','selF','selFD1','selFD2','selIdx','lat','disGating','lw','selD1','selD2','selD1D2','lin'};
    statNames={'All','Target Selective (%)','Distractor Selective (%)','D1 Selective (%)','D2 Selective (%)','CS','LMS','NMS (%)',...
        'Receptive field size','Memory field size (D1)','Memory field size (D2)','Mixed selectivity',...
        'Selectivity strength','Selectivity strength (D1)','Selectivity strength (D2)','Stimulus selectivity index','Response latency (s)',...
        'Distractor filtering (%)','Memory/Motor ratio','D1 selectivity index','D2 selectivity index','D1 & D2 selectivity index','Dummy linear'};
    %statNames={'Mean','D1resp','D2resp','Sel','D2sel','CS','LMS','NMS',...
    %    'D1D2rf','D1rf','D2rf','Mixedstr','D1D2pev','D1pev','D2pev','selIdx','lat','disGating','lw','selD1','selD2','selD1D2'};
    regions={'sDLPFC','iDLPFC'};




    for regionIdx=RegionsUsed
        region=regions{regionIdx};
        switch region
            case {'sDLPFC'}
                %[a8AD,a946D,Pbreak1AllD,Pbreak2AllD,oneBreakCont,~]=extract2Breaks(sDLPFCMeasuresCont,sDLPFCstatsUsed,statNames,region,MonkyID(end-1));
                [a8AD,a946D,Pbreak1AllD,Pbreak2AllD]=extract2Breaks(sDLPFCMeasuresCont,sDLPFCstatsUsed,statNames,region,MonkyID(end-1));

                % sDLPFC plot
                figure
                subplot(1,9,[2:4]);hold on

                %lims
                yticks([0:numel(order)]);set(gca,'yticklabel',statNames([1,order]));
                ytickangle(10);
                addSkippedTicks([-4:2:14],'x')
                xlim([-4 14.001])
                ylim([0 0+numel(order)])
                ylims=ylim;

                % draw boundaries and jitter box
                area([2.8-jitterBoundary 2.8+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                area([9.5-jitterBoundary 9.5+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                line([2.8 2.8],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)
                line([9.5 9.5],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)

                % xyt labels
                xlabel('AP position (mm)','FontWeight','bold')
                ylabel('Functional measures','FontWeight','bold')
                title('dDLPFC (Monkey 2)','FontWeight','normal')
                yBreaks=[Pbreak1AllD(:,2);Pbreak2AllD(:,2)];
                xMeasure_PD=[Pbreak1AllD(:,1);Pbreak2AllD(:,1)];
                [~,yAxOrder] = ismember(yBreaks,order);
                scatter(xMeasure_PD,yAxOrder,80,[0 0 0],'s','filled','MarkerEdgeColor','k','LineWidth',2)
                scatter([a8AD; a946D],[0 0],80,'r','s','filled','MarkerEdgeColor','k','LineWidth',2)
                upFontSize(18,0.012)
                ylim([0 0+numel(order)])
                saveFigure('em/gradient/breakPositions/',['SurrV2_' MonkyID 'breakplot'],'');
                M1break1AllD=Pbreak1AllD(:,1);
                save('~/Desktop/HPC_PK/data/em/gradient/SurrV2_M1breaks.mat','M1break1AllD');
            case {'iDLPFC'}

        end
    end
else
    MonkyID='M1_';
    RegionsUsed=[1 2]; %all s and iDLPFC
    sDLPFCstatsUsed=order;
    iDLPFCstatsUsed=order;
    %% Load data and set params
    % Load single monkey data
    load(['~/Desktop/HPC_PK/data/em/3regions_' MonkyID 'breakpoints.mat'])
    %load(['~/Desktop/HPC_PK/data/em/SurrV2_3regions_' MonkyID '_rbreakpointMeasuresCont_p005.mat'])
    %load(['~/Desktop/HPC_PK/data/em/3Bregion_surr_' MonkyID '_' statType 'breakpointMeasuresCont_p005.mat'])
    %
    oneBreakCont=[];
    twoBreakCont=[];
    a8AD=NaN;
    a946D=NaN;
    a8AV=NaN;
    a946V=NaN;
    statNames={'Mean','D1resp','D2resp','Sel','D2sel','CS','LMS','NMS',...
        'D1D2rf','D1rf','D2rf','Mixedstr','D1D2pev','D1pev','D2pev','selIdx','lat','disGating','lw','selD1','selD2','selD1D2'};
    statNames={'All','Target Selective (%)','Distractor Selective (%)','D1 Selective (%)','D2 Selective (%)','CS','LMS','NMS (%)',...
        'Receptive field size','Memory field size (D1)','Memory field size (D2)','Mixed selectivity',...
        'Selectivity strength','Selectivity strength (D1)','Selectivity strength (D2)','Stimulus selectivity index','Response latency (s)',...
        'Distractor filtering (%)','Memory/Motor ratio','D1 selectivity index','D2 selectivity index','D1 & D2 selectivity index','Dummy linear'};
    regions={'sDLPFC','iDLPFC'};
    for regionIdx=RegionsUsed
        region=regions{regionIdx};
        switch region
            case {'sDLPFC'}
                [a8AD,a946D,Jbreak1AllD,Jbreak2AllD]=extract2Breaks(sDLPFCMeasuresCont,sDLPFCstatsUsed,statNames,region,MonkyID(end-1));

                % sDLPFC plot
                figure
                subplot(1,9,[2:4]);hold on

                %lims
                yticks([0:numel(order)]);set(gca,'yticklabel',statNames([1,order]));
                ytickangle(10);
                addSkippedTicks([-4:2:14],'x')
                xlim([-4 14.001])
                ylim([0 0+numel(order)])
                ylims=ylim;

                % draw boundaries and jitter box
                area([2.8-jitterBoundary 2.8+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                area([9.5-jitterBoundary 9.5+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                line([2.8 2.8],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)
                line([9.5 9.5],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)

                % xyt labels
                xlabel('AP position (mm)','FontWeight','bold')
                ylabel('Functional measures','FontWeight','bold')
                title('dDLPFC (Monkey 1)','FontWeight','normal')
                yBreaks=[Jbreak1AllD(:,2);Jbreak2AllD(:,2)];
                xMeasure_JD=[Jbreak1AllD(:,1);Jbreak2AllD(:,1)];
                %[~,]=find(yBreaks==order);
                [~,yAxOrder] = ismember(yBreaks,order);
                scatter(xMeasure_JD,yAxOrder,80,[0 0 0],'s','filled','MarkerEdgeColor','k','LineWidth',2)
                scatter([a8AD; a946D],[0 0],80,'r','s','filled','MarkerEdgeColor','k','LineWidth',2)
                upFontSize(18,0.012)
                ylim([0 0+numel(order)])

            case {'iDLPFC'}
                [a8AV,a946V,Jbreak1AllV,Jbreak2AllV]=extract2Breaks(iDLPFCMeasuresCont,iDLPFCstatsUsed,statNames,region,MonkyID(end-1));

                % sDLPFC plot
                subplot(1,9,[6:8]);hold on

                %lims
                yticks([0:numel(order)]);set(gca,'yticklabel',[]);
                ytickangle(10);
                addSkippedTicks([-4:2:14],'x')
                xlim([-4 14.001])
                ylim([0 0+numel(order)])
                ylims=ylim;

                % draw boundaries and jitter box
                area([2.8-jitterBoundary 2.8+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                area([9.5-jitterBoundary 9.5+jitterBoundary],[ylims(2) ylims(2)], 'LineStyle',':', 'FaceAlpha',0.3, 'EdgeAlpha',0, 'FaceColor',[.85 .85 .85]);
                line([2.8 2.8],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)
                line([9.5 9.5],ylim,'color',[0.75 0.75 0.75],'lineStyle','--','lineWidth',2.5)

                % xyt labels
                %xlabel('AP position (mm)','FontWeight','bold')
                %ylabel('Functional measures','FontWeight','bold')
                title('vDLPFC (Monkey 1)','FontWeight','normal')
                yBreaks=[Jbreak1AllV(:,2);Jbreak2AllV(:,2)];
                xMeasure_JV=[Jbreak1AllV(:,1);Jbreak2AllV(:,1)];
                [~,yAxOrder] = ismember(yBreaks,order);
                scatter(xMeasure_JV,yAxOrder,80,[0 0 0],'s','filled','MarkerEdgeColor','k','LineWidth',2)
                scatter([a8AV; a946V],[0 0],80,'r','s','filled','MarkerEdgeColor','k','LineWidth',2)
                upFontSize(18,0.012)

                ylim([0 0+numel(order)])

                saveFigure('em/gradient/breakPositions/',['SurrV2_' MonkyID 'breakplot'],'');
                M2break1AllD=Jbreak1AllD(:,1);
                M2break2AllD=Jbreak2AllD(:,1);
                M2break1AllV=Jbreak1AllV(:,1);
                save('~/Desktop/HPC_PK/data/em/gradient/SurrV2_M2breaks.mat','M2break1AllD','M2break2AllD','M2break1AllV');
        end
    end
end

%[prctile8Av,prctile8Ad,prctile946d]=getFuncBoundPrctile(xMeasure_PD,Jbreak1AllD(:,1),Jbreak2AllD(:,1));

end