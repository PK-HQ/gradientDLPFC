function [yStr,yLimitsdDLFPC,yTickIntervaldDLPFC,yLimitsvDLFPC,yTickIntervalvDLPFC,funcMeasureStr,yLabels,xLabel]=getYplotparam(Monky,combineAx)

funcMeasureStr={'All','Tarselpct','Disselpct','D1selpct','D2selpct','CS','LMS','NMSpct',...
    'Trf','D1rf','D2rf','Mixedstr','selF','selFD1','selFD2','selIdx','lat','disGating','lw','selD1','selD2','selD1D2','lin'};
yLabels={'All','Target Selective (%)','Distractor Selective (%)','D1 Selective (%)','D2 Selective (%)','CS','LMS','NMS (%)',...
    'Receptive field size','Memory field size (D1)','Memory field size (D2)','Mixed selectivity',...
    'Selectivity strength','Selectivity strength (D1)','Selectivity strength (D2)','Stimulus selectivity index','Response latency (s)',...
    'Distractor filtering (%)','Memory/Motor ratio','D1 selectivity index','D2 selectivity index','D1 & D2 selectivity index','Dummy linear'};
xLabel={'AP location (mm)'};

%P
ylimContSup1=[0 99;50 150;50 125;0 125;0 125;0 125;0 125;0 125;0 8;0 8;0 8;0 8;0 50;0 50;0 25;0 .5;0 1.0;100 200;0.5 2;0 .5;
    0 .6;0 .6;0 14];

ylimSpacingContSup1=[33;25;25;25;25;25;25;25;1;1;1;1;10;5;5;0.1;.1;50;0.25;.1;.1;.1;2];

ylimContInf1=[0 99;50 150;50 125;0 125;0 125;0 125;0 125;0 125;0 6;0 6;0 6;0 8;0 50;0 30;0 30;0 .5;0 0.8;100 225;0.5 2;0 .5;
    0 .4;0 .5;0 14];

ylimSpacingContInf1=[33;25;25;25;25;25;25;25;1;1;1;1;10;5;5;0.1;.1;50;0.25;.1;.1;.1;2];

%J
ylimContSup2=[0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 7;0 7;0 7;0 2;0 50;0 10;0 7;0 .4;0 1.8;50 450;0.5 2;0 .4;0 .4;0 .4;0 14];

ylimSpacingContSup2=[0;0;0;0;0;0;0;0;1;1;1;.5;1;2;1;0.05;.2;50;0.25;.1;.1;.1;2];

ylimContInf2=[0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 8;0 8;0 8;0 8;0 50;0 25;0 18;0 .4;0 1.0;100 500;0.5 2;0 .4;0 .4;0 .4;0 14];

ylimSpacingContInf2=[0;0;0;0;0;0;0;0;1;1;1;2;1;4;3;0.05;.1;50;0.25;.1;.1;.1;2];

if combineAx==0
    if Monky==1 %P
        yStr={' iDLPFC (M2)                  sDLPFC (M2)'};
        yLimitsdDLFPC=ylimContSup1;

        yTickIntervaldDLPFC=ylimSpacingContSup1;
    
        yLimitsvDLFPC=ylimContInf1;

        yTickIntervalvDLPFC=ylimSpacingContInf1;

    elseif Monky==2 %J
        yStr={' iDLPFC (M1)                  sDLPFC (M1)'};
        yLimitsdDLFPC=ylimContSup2;

        yTickIntervaldDLPFC=ylimSpacingContSup2;
    
        yLimitsvDLFPC=ylimContInf2;

        yTickIntervalvDLPFC=ylimSpacingContInf2;
    end
elseif combineAx==1
    if Monky==1
        yStr={' vDLPFC (M2)                  dDLPFC (M2)'};
        yLimitsdDLFPC=[min(ylimContSup1(:,1),ylimContSup2(:,1)),max(ylimContSup1(:,2),ylimContSup2(:,2))];
        yTickIntervaldDLPFC=max(ylimSpacingContSup1(:,1),ylimSpacingContSup2(:,1));
        yLimitsvDLFPC=[min(ylimContInf1(:,1),ylimContInf2(:,1)),max(ylimContInf1(:,2),ylimContInf2(:,2))];
        yTickIntervalvDLPFC=max(ylimSpacingContInf1(:,1),ylimSpacingContInf2(:,1));
    elseif Monky==2
        yStr={' vDLPFC (M1)                  dDLPFC (M1)'};
        yLimitsdDLFPC=[min(ylimContSup1(:,1),ylimContSup2(:,1)),max(ylimContSup1(:,2),ylimContSup2(:,2))];
        yTickIntervaldDLPFC=max(ylimSpacingContSup1(:,1),ylimSpacingContSup2(:,1));
        yLimitsvDLFPC=[min(ylimContInf1(:,1),ylimContInf2(:,1)),max(ylimContInf1(:,2),ylimContInf2(:,2))];
        yTickIntervalvDLPFC=max(ylimSpacingContInf1(:,1),ylimSpacingContInf2(:,1));
    end
end

