function addSkippedTicks(tickLabels,axisType)
switch class(tickLabels(1))
    case 'double'
        switch axisType
            case {'x'}
                toSkipIdx=2:2:numel(tickLabels);
                tickLabelCells=num2cell(tickLabels);
                tickLabelCells(toSkipIdx)={' '};
                set(gca,'xtick',tickLabels,'xticklabel',tickLabelCells);    
            case {'y'}
                toSkipIdx=2:2:numel(tickLabels);
                tickLabelCells=num2cell(tickLabels);
                tickLabelCells(toSkipIdx)={' '};
                set(gca,'ytick',tickLabels,'yticklabel',tickLabelCells);
        end
    case 'cell'
        names = ({'Tar'; ''; 'Dis'; '';'Cue'});
        set(gca,'xtick',1:5,'xticklabel',names);
end