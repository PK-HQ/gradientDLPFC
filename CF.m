function CF()
%% Closes all plotted figures
    % find all handles of axes (graphs) 
    axh = findall(groot,'type','axes');
    % get handles of parent figures containing graphs
    fxh = get(axh,'parent');
    switch isa(fxh,'cell')
        case {1}
            %// close figures containg axes
            close(fxh{:});
        case {0}
            %
    end
end