function [comp1,comp2,a,b,minmi,mi] = Decorrelation(delay1,delay2)
% Input: both delay1 and delay2 are of size N x loc, where N is the number of
% neuron, and loc is the size of target locations. They contain the trial
% and time averaged data (second half of each delay) in Delay 1 and Delay 2

% Output: comp1 and comp2 are of size N x loc, and has minimal mutual
% inforamtion. 
% minmi: the mutual information
% a and b: the mixing coefficients


arange = 0:0.01:0.99; % range for parameter search
brange = 0:0.01:0.99;

mi = zeros(length(arange),length(brange)); % mutual information
for aind = 1:length(arange)
    for bind = 1:length(brange)
        a = arange(aind);
        b = brange(bind);
        
        x = (b*delay1(:)-delay2(:))/(a*b-1);% preparation component
        y = (a*delay2(:)-delay1(:))/(a*b-1);% memory component
        mi(aind,bind) = computeMI(x,y);
    end
end
[I,J] = find(mi == min(min(mi)));
a = arange(I);
b = brange(J);
minmi = min(min(mi));

comp1 = (delay1-a*delay2)/(1-a*b);
comp2 = (delay2-b*delay1)/(1-a*b);
end


function mi = computeMI(x,y)
binCount = 7;
xgrid = zeros(binCount+1,1);
ygrid = zeros(binCount+1,1);
mi = 0;
binWidth = 1/binCount*100;

for i = 1:binCount+1
    xgrid(i) = prctile(x,min(100,(i-1)*binWidth));
    ygrid(i) = prctile(y,min(100,(i-1)*binWidth));
end


label = zeros(length(x),2);
for i = 1:length(x)
    % for x label
    for j = 1:binCount
        if xgrid(j)<=x(i)&&x(i)<xgrid(j+1)
            label(i,1) = j;
            break;
        end
    end
    % for y label
    for j = 1:binCount
        if ygrid(j)<=y(i)&&y(i)<ygrid(j+1)
            label(i,2) = j;
            break;
        end
    end
end

ss = size(label,1); % sample size

for i = 1:binCount
    for j = 1:binCount
        P_ij = length(find(label(:,1)==i&label(:,2)==j))/ss;
        if P_ij~=0
            P_i = length(find(label(:,1)==i))/ss;
            P_j = length(find(label(:,2)==j))/ss;
            mi = mi + P_ij*log(P_ij/(P_i*P_j));
        end
    end
end


end


