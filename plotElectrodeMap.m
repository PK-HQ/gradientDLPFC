function plotElectrodeMap(folderPathEMPlot,referenceElectrodeMap,...
    desiredStatistic)
electrodeSize=27;

%end boundaries for each region
a6DR=3.5;
a8B=7.9;
a9=13.9;
a10SUP=14.5;
FPsup=16.4;

a8AD=2.8;
a946D=9.5;
a46D=13.5;
a10D=15.9;
FPd=21.3;

a8AV=3;
a946V=9.7;
a46V=13.8;
a10V=15.9;
FPv=21.3;

%% Reference electrode map 632 x 6
%   monkey |  arrayNo |     arrayAB    |  electrodeNo  |     unitNo      |  sessionNo
%    P=1   |   1-12   |  A=1 B=2 nil=0 |      1-32     |  1-5 (useless)  |     1-8
%anatLocations={'vFEF';'dFEF';'v46';'dDLPFC';'8b';'SEF';'vDLPFC';'d46';'a8';'fef';'dDLPFC';'vDLPFC'};
%                1        2     3       4      5    6       7      8    9     10     11        12 
xLab={'AP location (mm)'};

%                1        2     3       4      5    6       7      8    9     10    11     12       13
%fprintf('Plotting %s\n',statNames{stat})
%merge each electrode with its statistic of interest
%monkeysDataRef=[referenceElectrodeMap,refAllNeurons];
monkeysData=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,7:8)];


%region-array parcels
dDLPFC=[1 2 4 8 10 11 12];%[2 4 8]; %[2 4 8 10 11 12];
vDLPFC=[3 7];%[1 3 7]; %[1 3 7 13];
PMA=[5 6 9];%;[5 6 9]; %[5 6 9];

dDLPFC=[1 2 4 8 10 11 12];%[2 4 8]; %[2 4 8 10 11 12];
vDLPFC=[3 7 13];%[1 3 7]; %[1 3 7 13];
PMA=[5 6 9];%;[5 6 9]; %[5 6 9];

wantedregions={'premotor area','dDLPFC','vDLPFC','LPFC'};

[greenOut,magentaOut]=getColors(9,1);
colors=[greenOut(4,:); greenOut(6,:); greenOut(8,:)];

%% plotting
for monky=[1 2]
    switch monky
        case {1}
            monkyStr='2';
            a8AD=2.6222;
            a946D=NaN;
            a46D=NaN;
            a10D=NaN;
            a8AV=2.6222*(3.4333/3.4556);
            a946V=NaN;
            a46V=NaN;
            a10V=NaN;
        case {2}
            monkyStr='1';
            a8AD=3.4556;
            a946D=9.2;
            a46D=NaN;
            a10D=NaN;
            a8AV=3.4333;
            a946V=9.2;
            a46V=NaN;
            a10V=NaN;
    end
    for reg=4
        wantedregion=wantedregions{reg};fprintf('%s\n',wantedregion);

        %% define wanted region and breakpoints
        switch wantedregion
            case {'premotor area'}
                [regionNeurons,~,~]=find(monkeysData(:,2)==PMA);

                breakpoints=[a6DR,a8B,a9,a10SUP];
                %define breakpoints
                v1=breakpoints(1);
                v2=breakpoints(2);
                v3=breakpoints(3);
                v4=breakpoints(4);

                dataColor=colors(1,:);

            case {'dDLPFC'}
                [regionNeurons,~,~]=find(monkeysData(:,2)==dDLPFC);

                breakpoints=[a8AD,a946D,a46D,a10D];
                %define breakpoints
                v1=breakpoints(1);
                v2=breakpoints(2);
                v3=breakpoints(3);
                v4=breakpoints(4);
                dataColor=colors(2,:);


            case {'vDLPFC'}
                [regionNeurons,~,~]=find(monkeysData(:,2)==vDLPFC);


                breakpoints=[a8AV,a946V,a46V,a10V];
                %define breakpoints
                v1=breakpoints(1);
                v2=breakpoints(2);
                v3=breakpoints(3);
                v4=breakpoints(4);
                dataColor=colors(3,:);

            case {'LPFC'}
                [regionNeurons,~,~]=find(monkeysData(:,2)==[PMA dDLPFC vDLPFC]);

        end
        load('~/Desktop/HPC_PK/data/neuronidx/allElectrodeCoords.mat')
        monkeysData1=[zeros(size(allElectrodeCoords,1),5),allElectrodeCoords];
        monkeysData1=monkeysData1(monkeysData1(:,end-3)==monky,:);
        %% Define X and Y
        %monkeysData1=monkeysData(regionNeurons,:);
        x=[monkeysData1(:,end-1)];
        y=[monkeysData1(:,end)];
        desiredStatistic=zeros(size(monkeysData1(:,7),1)); %rand(size(monkeysData1(:,7),1),1);
        regionsId=[x monkeysData1(:,2)];

        %get mean stats of each unique electrode positions
        meanUniqPosDatas=collapsePosition([x y desiredStatistic],'mean');
        x1=meanUniqPosDatas(:,1);
        x2=meanUniqPosDatas(:,2);
        y=meanUniqPosDatas(:,3);



    %% 2d gscatter
        % nan-zscore
        [y,~,~]=nanzscore(y);
        %rm nan
        x1=x1(~isnan(y));
        x2=x2(~isnan(y));
        y=y(~isnan(y));

        %def x, y, stat
        data=roundData(y,.5);
        range=(numel(unique(data)));

        %plot scatter
        figure
        colorMap=colormap(cubehelix(256,0.43,0.35,2.5,1,[0,.95],[0,1]));
        colorMap=colorMap(1:26:256,:);hold on
        gscatter(x1,x2,zeros(size(data)),[0 0 0; colorMap(3:end,:)],'.',electrodeSize+12);
        gscatter(x1,x2,data,[.7 .7 .7; colorMap(3:end,:)],'.',electrodeSize);

        %labels
        ylabel('DV location (mm)');
        xlabel(xLab{1});
        title(['Monkey ' monkyStr])
        set(gca,'Color','w')
        axis equal
        %hleg=get(gca,'Legend');
        %title(hleg, 'Z-scores');


        % region bounds
        xlim([-5 15]);xticks(-6:2:16);
        ylim([-5 11]);yticks(-6:2:12);

        % anat boundaries
        line([a6DR a6DR],[8 15],'color','k','lineStyle','--','lineWidth',2)
        line([a8B a8B],[8 15],'color','k','lineStyle','--','lineWidth',2)
        line([a9 a9],[8 15],'color','k','lineStyle','--','lineWidth',2)

        line([a8AD a8AD],[0 8],'color','k','lineStyle','--','lineWidth',2)
        line([a946D a946D],[0 8],'color','k','lineStyle','--','lineWidth',2)
        %line([a46D a46D],[0 8],'color','k','lineStyle','--','lineWidth',2)

        line([a8AV a8AV],[-6 0],'color','k','lineStyle','--','lineWidth',2)
        line([a946V a946V],[-6 0],'color','k','lineStyle','--','lineWidth',2)
        %line([a46V a46V],[-6 0],'color','k','lineStyle','--','lineWidth',2)


        line([-6 25],[0 0],'color','k','lineStyle','--','lineWidth',2)
        upFontSize(20,0.02)
        
        hold off;


        %% for poster
        ylim([-7 8]);xlim([-5 14])

    end
    [hleg,~]=legend('show');
    title(hleg,'Z')
    upFontSize(30,0.013)
    set(gca,'linewidth',3.5)
    removeURTicks
    saveFigure(folderPathEMPlot,['monkey_' monkyStr '_anat'],'')

end

    
end
        
        