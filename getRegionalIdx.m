function getRegionalIdx(folderPathNeuronIdx,referenceElectrodeMap,desiredStatistic,anatLocations)
versionStr='BothMonkeys';
%% old 
a6DR=3.5;
a8B=7.9;

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
monkyStrs={'M1','M2'};
%% new
%{
a8AD=2.70;
a946D=9.5;
a46D=13.5;
a10D=15.9;
FPd=21.3;

a8AV=2.2;
a946V=NaN;
a46V=13.8;
a10V=15.9;
FPv=21.3;
%}

for Monky=[1 2]
    monkyStr=monkyStrs{Monky};
    %[a8AD,a946D,a8AV,a946V]=getFunctionalBounds()
    [a8AD,a946D,a8AV,a946V]=getFunctionalBoundsIndiv(Monky);
    if Monky==1
        a946D=12;
    end
    %a8AD=2.4;
    %a946D=20;
    %a8AV=99;
    %a946V=99;

    %a8AD=3.42;
    %a946D=9.5;
    %a8AV=3.43;
    %a946V=9.2;
    
    %% Reference electrode map 632 x 6
    %   monkey |  arrayNo |     arrayAB    |  electrodeNo  |     unitNo      |  sessionNo
    %    P=1   |   1-12   |  A=1 B=2 nil=0 |      1-32     |  1-5 (useless)  |     1-8
    %anatLocations={'vFEF';'dFEF';'v46';'dDLPFC';'8b';'SEF';'vDLPFC';'d46';'a8';'fef';'dDLPFC';'vDLPFC'};
    %                1        2     3       4      5    6       7      8    9     10     11        12 
    monkeysData=[referenceElectrodeMap(:,1:6) desiredStatistic referenceElectrodeMap(:,7:8)];
    yCol=9;
    xCol=8;

    %get neurons in each region
    [regionNeuronsa6DR,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)>8 & monkeysData(:,xCol)<=a6DR); %ylims
    [regionNeuronsa8B,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)>8 & monkeysData(:,xCol)>a6DR); %ylims

    [regionNeuronsa8AD,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)>0 & monkeysData(:,yCol)<=8 & monkeysData(:,xCol)<=a8AD); %ylims
    [regionNeuronsa946D,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)>0 & monkeysData(:,yCol)<=8 & monkeysData(:,xCol)>a8AD & monkeysData(:,xCol)<=a946D); %ylims
    [regionNeuronsa46D,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)>0 & monkeysData(:,yCol)<=8 & monkeysData(:,xCol)>a946D); %ylims

    [regionNeuronsa8AV,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)<=0 & monkeysData(:,xCol)<=a8AV); %ylims
    [regionNeuronsa946V,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)<=0 & monkeysData(:,xCol)>a8AV); %ylims
    [regionNeuronsa46V,~,~]=find(monkeysData(:,1)==Monky & monkeysData(:,yCol)<=0 & monkeysData(:,xCol)>a946V); %ylims

    %put into regionalNeuronsIdx
	regionalNeuronsIdx=zeros(632,8);
    regionalNeuronsIdx(1:numel(regionNeuronsa8AV),1)...
        =regionNeuronsa8AV;
    regionalNeuronsIdx(1:numel(regionNeuronsa8AD),2)...
        =regionNeuronsa8AD;
    regionalNeuronsIdx(1:numel(regionNeuronsa46V),3)...
        =regionNeuronsa46V;
    regionalNeuronsIdx(1:numel(regionNeuronsa946D),4)...
        =regionNeuronsa946D;
    regionalNeuronsIdx(1:numel(regionNeuronsa8B),5)...
        =regionNeuronsa8B;
    regionalNeuronsIdx(1:numel(regionNeuronsa6DR),6)...
        =regionNeuronsa6DR;
    regionalNeuronsIdx(1:numel(regionNeuronsa946V),7)...
        =regionNeuronsa946V;
    regionalNeuronsIdx(1:numel(regionNeuronsa46D),8)...
        =regionNeuronsa46D;

    %print change in n-neuron from previous
    fprintf('vFEF %.0f\ndFEF %.0f\nv46 %.0f\ndDLPFC %.0f\n8b %.0f\n6DR %.0f\nvDLPFC %.0f\nd46 %.0f\n',...
        numel(regionNeuronsa8AV),numel(regionNeuronsa8AD),...
        NaN,numel(regionNeuronsa946D),...
        numel(regionNeuronsa8B),numel(regionNeuronsa6DR),...
        numel(regionNeuronsa946V),numel(regionNeuronsa46D));
    numel(nonzeros(regionalNeuronsIdx))
    %save regionalNeuronsIdx
    save([folderPathNeuronIdx monkyStr '_regionalNeuronsIdx_' versionStr '.mat'],'regionalNeuronsIdx')

end
%remove duplicate cells, combine M1 and M2 cells by regions
removeDuplicatesCombineMonkeyCells(folderPathNeuronIdx,versionStr,anatLocations)

end





