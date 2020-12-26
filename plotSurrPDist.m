function plotSurrPDist(sDLPFCsurrPMeasuresCont,iDLPFCsurrPMeasuresCont,Monky)
%% TESTS IF DATA DISTRIBUTIONS HAVE NORMAL OR WEIBULL FIT, PLOTS ALL TESTED FITS AND BEST FITS

titlestr='all periods';
str='all';
p1Normcont=[];
p1Weibullcont=[];
p2Normcont=[];
p2Weibullcont=[];

%colors
B=[115,115,115;
82,82,82;
37,37,37;
158,202,225;
66,146,198;
8,81,156];

P=[158,154,200;
117,107,177;
84,39,143];

R=[251,106,74;
222,45,38;
165,15,21];

O=[115,115,115;
82,82,82;
37,37,37;
253 165 15;
248 222 126;
255 211 0;];
textR=[165,15,21]*1.3/255; textB=[8,81,220]*0.5/255;
%% Process raw
%sDLPFC
if Monky==1
    MonkyStr='2';
elseif Monky==2
    MonkyStr='1';
elseif Monky==3
    MonkyStr='1&2';
end
pvalues1seg=nonan(cell2mat(sDLPFCsurrPMeasuresCont(:,1)));
pvalues2seg=nonan(cell2mat(sDLPFCsurrPMeasuresCont(:,2)));
pvalues3seg=nonan(cell2mat(sDLPFCsurrPMeasuresCont(:,3)));

figure
subplot(7,1,[1 2 3]);
%% 1seg
%Get weibull fit parameters
meanDefault=0.5;
mean1seg=mean(pvalues1seg);
mean2seg=mean(pvalues2seg);

weibullParams=[0.5 0.25; 0.5 0.5; 0.5 0.75]; %exponential with different steepness
normParams1=[meanDefault 0.1; meanDefault 0.15; meanDefault 0.30; mean1seg 0.1; mean1seg 0.15; mean1seg 0.30]; %exponential with different steepness
normParams2=[meanDefault 0.1; meanDefault 0.15; meanDefault 0.30; mean2seg 0.1; mean2seg 0.15; mean2seg 0.30];
%% STATS: 1seg weibull/normal test
%normal
for normParamSet=1:length(normParams1)
    normParam=normParams1(normParamSet,:);
    [~,p1,adstat1]=adtestCustom(pvalues1seg,normParam,[],'Distribution','norm');
    if p1>.05
        sprintf('1-seg deg data may be norm shape %.2f, %.2f (p=%.2f)',normParam(1),normParam(2),p1)
    elseif p1<.05
        %sprintf('Data NOT norm shape %.2f, %.2f (p=%.2f)',normParam(1),normParam(2),p1)
    end
    %store in cont
    p1Normcont(normParamSet,1)=p1;
    p1Normcont(normParamSet,2)=adstat1;
    p1Normcont(normParamSet,3)=normParam(1);
    p1Normcont(normParamSet,4)=normParam(2);
end
%weibull
for weibullParamSet=1:length(weibullParams)
    weibullParam=weibullParams(weibullParamSet,:);
    [~,p1,adstat1]=adtestCustom(pvalues1seg,[],weibullParam,'Distribution','weibull'); 
    if p1>.05
        sprintf('1-seg data may be weibull shape %.1f, %.1f (p=%.2f)',weibullParam(1),weibullParam(2),p1)
    elseif p1<.05
        %sprintf('Data NOT weibull shape %.1f, %.1f (p=%.2f)',weibullParam(1),weibullParam(2),p1)
    end
    %store in cont
    p1Weibullcont(weibullParamSet,1)=p1;
    p1Weibullcont(weibullParamSet,2)=adstat1;
    p1Weibullcont(weibullParamSet,3)=weibullParam(1);
    p1Weibullcont(weibullParamSet,4)=weibullParam(2);
end

%% STATS: 2seg weibull/normal test
%normal
for normParamSet=1:length(normParams1)
    normParam=normParams2(normParamSet,:);
    [~,p2,adstat2]=adtestCustom(pvalues2seg,normParam,[],'Distribution','norm');
    if p2>.05
        sprintf('2-seg data may be norm shape %.2f, %.2f (p=%.2f)',normParam(1),normParam(2),p2)
    elseif p2<.05
        %sprintf('Data NOT norm shape %.2f, %.2f (p=%.2f)',normParam(1),normParam(2),p2)
    end
    %store in cont
    p2Normcont(normParamSet,1)=p2;
    p2Normcont(normParamSet,2)=adstat2;
    p2Normcont(normParamSet,3)=normParam(1);
    p2Normcont(normParamSet,4)=normParam(2);
end
%weibull
for weibullParamSet=1:length(weibullParams)
    weibullParam=weibullParams(weibullParamSet,:);
    [~,p2,adstat2]=adtestCustom(pvalues2seg,[],weibullParam,'Distribution','weibull');
    if p2>.05
        sprintf('2-seg data may be weibull shape %.1f, %.1f (p=%.2f)',weibullParam(1),weibullParam(2),p2)
    elseif p2<.05
        %sprintf('Data NOT weibull shape %.1f, %.1f (p=%.2f)',weibullParam(1),weibullParam(2),p2)
    end
    %store in cont
    p2Weibullcont(weibullParamSet,1)=p2;
    p2Weibullcont(weibullParamSet,2)=adstat2;
    p2Weibullcont(weibullParamSet,3)=weibullParam(1);
    p2Weibullcont(weibullParamSet,4)=weibullParam(2);
end
header={'P-val','AD-stat','Param 1','Param 2'};
p1NormcontH=[header;num2cell(p1Normcont(3:6,:))];
p2NormcontH=[header;num2cell(p2Normcont(3:6,:))];
p1WeibullcontH=[header;num2cell(p1Weibullcont)];
p2WeibullcontH=[header;num2cell(p2Weibullcont)];

%% PLOT: All tested weibull/norm distributions, for 1-seg and 2-seg
figure; subplot(7,1,[1 2 3]); hold on;
for dist=1:size(p1Normcont,1)
    lineColor=B(dist,:)/255;
    normal1=makedist('Normal','mu',p1Normcont(dist,3),'sigma',p1Normcont(dist,4));
    xval=0:0.001:1;
    normY1=pdf(normal1,xval);
    plot(xval,normY1,'color',lineColor,'LineWidth',3);

end
for dist=1:size(p1Weibullcont,1)
    lineColor=P(dist,:)/255;
    weibull2=makedist('Weibull','a',p1Weibullcont(dist,3),'b',p1Weibullcont(dist,4));
    xval=0:0.001:1;
    weibullY2=pdf(weibull2,xval);
    plot(xval,weibullY2,'color',lineColor,'LineWidth',3);
end
hold off;
ylim([0 8])
xlim([0 1])
ylims=ylim;
legend({'?=0.50, ?=0.10','?=0.50, ?=0.15','?=0.50, ?=0.30',...
    '?=0.45, ?=0.10','?=0.45, ?=0.15','?=0.45, ?=0.30',...
    'k=0.50, ?=0.25', 'k=0.50, ?=0.50', 'k=0.50, ?=0.75'},'location','eastoutside')
ylabel('Count','FontWeight','bold')
xlabel('p-values','FontWeight','bold')
title('Normal and Weibull p-value distributions, 1-seg data')
text(-1.16,ylims(2)+0.1*(ylims(2)-ylims(1)),'A','FontSize',20,'Fontweight','bold');
hold off
upFontSize(20,.01)
removeURTicks()

subplot(7,1,[5 6 7]); hold on;
for dist=1:size(p2Normcont,1)
    lineColor=O(dist,:)/255;
    normal1=makedist('Normal','mu',p2Normcont(dist,3),'sigma',p2Normcont(dist,4));
    xval=0:0.001:1;
    normY1=pdf(normal1,xval);
    plot(xval,normY1,'color',lineColor,'LineWidth',3);
end
for dist=1:size(p2Weibullcont,1)
    lineColor=R(dist,:)/255;
    weibull2=makedist('Weibull','a',p2Weibullcont(dist,3),'b',p2Weibullcont(dist,4));
    xval=0:0.001:1;
    weibullY2=pdf(weibull2,xval);
    plot(xval,weibullY2,'color',lineColor,'LineWidth',3);
end
ylim([0 8])
xlim([0 1])
ylims=ylim;
legend({'?=0.50, ?=0.10','?=0.50, ?=0.15','?=0.50, ?=0.30',...
    '?=0.36, ?=0.10','?=0.36, ?=0.15','?=0.36, ?=0.30',...
    'k=0.50, ?=0.25', 'k=0.50, ?=0.50', 'k=0.50, ?=0.75'},'location','eastoutside')
ylabel('Count','FontWeight','bold')
xlabel('p-values','FontWeight','bold')
title('Normal and Weibull p-value distributions, 2-seg data')
text(-1.16,ylims(2)+0.1*(ylims(2)-ylims(1)),'B','FontSize',20,'Fontweight','bold');
hold off
upFontSize(20,.01)
removeURTicks()

%% PLOT: Actual tested weibull/norm distributions
figure;subplot(7,1,[1 2 3 4]);
wantedFitIdxNorm1=find(p1Normcont(:,1)>0.05);
wantedFitIdxWeibull2=find(p2Weibullcont(:,1)>0.05);
wantedFitIdxNorm2=find(max(p2Normcont(3:6,2)))+3; %find highest AD, add 3 to get the highest of last 3 entries
wantedFitIdxWeibull1=find(max(p1Weibullcont(:,2))); %find highest AD
%Norm dist, 1-seg 
normal1=makedist('Normal','mu',p1Normcont(wantedFitIdxNorm1,3),'sigma',p1Normcont(wantedFitIdxNorm1,4));
xval=0:0.001:1;
normY1=pdf(normal1,xval);

% Weibull dist, 2-seg
weibull1=makedist('Weibull','a',p1Weibullcont(wantedFitIdxWeibull2,3),'b',p1Weibullcont(wantedFitIdxWeibull2,4));
xval=0:0.001:1;
weibullY1=pdf(weibull1,xval);

% Weibull dist, 2-seg
weibull2=makedist('Weibull','a',p2Weibullcont(wantedFitIdxWeibull2,3),'b',p2Weibullcont(wantedFitIdxWeibull2,4));
xval=0:0.001:1;
weibullY2=pdf(weibull2,xval);

% Normal dist, 2-seg
normal2=makedist('Normal','mu',p2Normcont(wantedFitIdxNorm1,3),'sigma',p2Normcont(wantedFitIdxNorm1,4));
xval=0:0.001:1;
normY2=pdf(normal2,xval);

%Plot here
plot(xval,normY1,'color','black','LineWidth',8);hold on;
plot(xval,weibullY1,'color','black','LineWidth',8);
plot(xval+0.0075,weibullY2,'color','black','LineWidth',8);
plot(xval,normY2,'color','black','LineWidth',8);

p1=plot(xval,normY1,'color',textB*1.1,'LineWidth',4); 
p1=plot(xval,weibullY1,'color',textB*1.1,'LineWidth',4);
p1=plot(xval+0.0075,weibullY2,'color',textR*1.1,'LineWidth',4);
p1=plot(xval,normY2,'color',textR*1.1,'LineWidth',4);


%data histograms, 1- and 2-seg
histogram(pvalues1seg,'binwidth',0.05,'FaceColor',textB,...
    'FaceAlpha',0.95,'LineWidth',2);
histogram(pvalues2seg,'binwidth',0.05,'FaceColor',textR,...
    'FaceAlpha',0.9,'LineWidth',2);

%XY-lims
ylim([0 4])
xlim([0 1])
ylims=ylim;
% median line
line([median(pvalues1seg),median(pvalues1seg)],ylim,'LineStyle','--','color',textB,'LineWidth',4)
line([median(pvalues2seg),median(pvalues2seg)],ylim,'LineStyle','--','color',textR,'LineWidth',4)
% text
text(0.6,ylims(2)-0.1*(ylims(2)-ylims(1)),['p_{1-seg, normal dist.}='...
    sprintf('%.2f',p1Normcont(wantedFitIdxNorm1,1))],'FontSize',6,...
    'FontWeight','Normal');
text(0.6,ylims(2)-0.2*(ylims(2)-ylims(1)),['p_{1-seg, Weibull dist.}'...
    sprintf('<0.0005')],'FontSize',6,...
    'FontWeight','Normal');
text(0.08,ylims(2)-0.1*(ylims(2)-ylims(1)),['p_{2-seg, normal dist.}'...
    sprintf('<0.0005')],'FontSize',6,...
    'FontWeight','Normal');
text(0.08,ylims(2)-0.2*(ylims(2)-ylims(1)),['p_{2-seg, Weibull dist.}='...
    sprintf('%.2f',p2Weibullcont(wantedFitIdxWeibull2,1))],'FontSize',6,...
    'FontWeight','Normal');

%% legend title
f=get(gca,'children');
f=flipud(f);
legend(f([6 7 9 10]),{'1-seg normal/Weibull fit','2-seg normal/Weibull fit','1-seg distribution','2-seg distribution'},'location','eastoutside')
title('P-value distribution of 1-seg and 2-seg models (n=12)','Fontweight','Normal')
%text(-0.075,ylims(2)+0.1*(ylims(2)-ylims(1)),'A','FontSize',20,'Fontweight','bold');
ylabel('Count','FontWeight','normal')
xlabel('Model-comparison p-value (p_{MC})','FontWeight','normal')
upFontSize(20,.01)
removeURTicks()
saveFigure('em',['surrPdist_TestedFitDists'],'')

end