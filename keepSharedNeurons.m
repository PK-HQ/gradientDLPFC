function contTrimmed=keepSharedNeurons(contToTrim,usedMeasures)
%% FIND CONSENSUS NEURONS SHARED ACROSS GIVEN MEASURES
presCont=[];
conts9to22=contToTrim(usedMeasures);

%find nueorns found across all measures that are being considered
for i=1:length(conts9to22)
    pres=~isnan(conts9to22{i});
    presCont=[presCont,pres];
end
   
%get the index
presentNeurons=sum(presCont,2)==size(presCont,2);
sum(presentNeurons);

%NaN all non-shared neurons
for j=usedMeasures
    toTrim=contToTrim{j};
    toTrim(~presentNeurons)=NaN;
    contToTrim{j}=toTrim;
end
contTrimmed=contToTrim;
end
   
