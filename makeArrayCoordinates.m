function [FMA36,FMA18,HDFMA18]=makeArrayCoordinates()
%% MAKES TEMPLATE ARRAYS WITH MILIMETER COORDINATES FOR THE 3 TYPES OF ARRAYS
%% ALL MEASUREMENTS ARE WRT THESE ORIGINS: 18FMA (ELECTRODE 18), 18HDFMA (ELECTRODE 15), 36FMA (ELECTRODE 36) 
%% PLEASE REFER TO MICROPROBES WEBSITE FOR THE FMA SCHEMATIC AND SPECIFICATIONS

%% 36 FMA (origin at electrode 36)
%dims of the space between edge of arary and electrode
xProtude=0.2;
xRecess=0.4;
yTop=0.3;
yMid=0.346;
yEnd=0.4620;
xSpacer=0.4;

rowP=xProtude+[0 (1:9-1)*xSpacer];
rowR=xRecess+[0 (1:9-1)*xSpacer];
FMA36_row1=[fliplr(28:36');rowP;repmat(yTop+yMid*0,numel(rowP),1)']'; %each row of electrodes
FMA36_row2=[fliplr(19:27');rowR;repmat(yTop+yMid*1,numel(rowR),1)']';
FMA36_row3=[fliplr(10:18');rowP;repmat(yTop+yMid*2,numel(rowP),1)']';
FMA36_row4=[fliplr(1:9');rowR;repmat(yTop+yMid*3,numel(rowR),1)']';
FMA36=[FMA36_row1;FMA36_row2;FMA36_row3;FMA36_row4];

%% 18 FMA (origin at electrode 18)
xProtude=0.4;
xRecess=0.6;
yTop=0.3;
yMid=0.346;
yEnd=0.612;
xSpacer=0.4;

rowP=xProtude+[0 (1:4)*xSpacer];
rowR=xRecess+[0 (1:3)*xSpacer];
FMA18_row1=[fliplr(15:18');rowR;repmat(yTop+yMid*0,numel(rowR),1)']';
FMA18_row2=[fliplr(10:14');rowP;repmat(yTop+yMid*1,numel(rowP),1)']';
FMA18_row3=[fliplr(6:9');rowR;repmat(yTop+yMid*2,numel(rowR),1)']';
FMA18_row4=[fliplr(1:5');rowP;repmat(yTop+yMid*3,numel(rowP),1)']';
FMA18=[FMA18_row1;FMA18_row2;FMA18_row3;FMA18_row4];

%% 18 HD FMA (origin at electrode 15)
xProtude=0.4;
xRecess=0.525;
yTop=0.3;
yMid=0.2165;
yEnd=0.6505;
xSpacer=0.25;

rowP=xProtude+[0 (1:4)*xSpacer];
rowR=xRecess+[0 (1:3)*xSpacer];
HDFMA18_row1=[fliplr(15:18');rowR;repmat(yTop+yMid*0,numel(rowR),1)']';
HDFMA18_row2=[fliplr(10:14');rowP;repmat(yTop+yMid*1,numel(rowP),1)']';
HDFMA18_row3=[fliplr(6:9');rowR;repmat(yTop+yMid*2,numel(rowR),1)']';
HDFMA18_row4=[fliplr(1:5');rowP;repmat(yTop+yMid*3,numel(rowP),1)']';
HDFMA18=[HDFMA18_row1;HDFMA18_row2;HDFMA18_row3;HDFMA18_row4];

end