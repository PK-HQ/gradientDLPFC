function meanUniqPosDatas=collapsePosition(data,collapseStat)
%% get mean of data at each unique x,y position

meanUniqPosDatas=[];
%find unique pos
uniquePositions = unique(data(:,1:2),'rows');

%mean data of each unique pos, store in 'Datas'
for row=1:length(uniquePositions)
    uniquePos=uniquePositions(row,:);
    %wantedRows=find(data(:,1:2)==uniquePos);
    %wantedRowsIdx=wantedRows(1:end-length(wantedRows)/2);

    %[~,index_A,~] = intersect(data(:,1:2),uniquePos,'rows');
    binaryIdx=data(:,1)==uniquePos(1) & data(:,2)==uniquePos(2);
    idx=find(binaryIdx,100);
    uniquePosData=data(idx,3);
    switch collapseStat
        case {'mean'}
            meanUniqPosData=mean(uniquePosData);
        case {'median'}
            meanUniqPosData=median(uniquePosData);
    end
    meanUniqPosDatas(row,1:3)=[uniquePos meanUniqPosData];
end
