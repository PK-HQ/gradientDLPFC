function saveFigure(folderPath,fileName,optLabel)
%SAVES FIG
mainPlotPath='/Volumes/Users/PK/Desktop/HPC_PK/plots/';

if ~exist([mainPlotPath folderPath '/eps/'], 'dir')
   mkdir([mainPlotPath folderPath '/eps/'])
end
if ~isempty(optLabel)
    fig = gcf;print(fig,[mainPlotPath folderPath '/' fileName '_' optLabel '.png'],'-dpng','-r450');
    fig = gcf;print(fig,[mainPlotPath folderPath '/eps/' fileName '_' optLabel '.eps'],'-depsc');
else 
    fig = gcf;print(fig,[mainPlotPath folderPath '/' fileName '.png'],'-dpng','-r450');
    fig = gcf;print(fig,[mainPlotPath folderPath '/eps/' fileName '.eps'],'-depsc');
end

end
