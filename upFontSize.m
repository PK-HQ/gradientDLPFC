function upFontSize(size,tickLength)
%INCREASE OVERALL FONT SIZE
    set(findall(gcf,'-property','FontSize'),'FontSize',size)
    set(gca,'linewidth',2.5)
    set(gca,'TickDir','out');
    set(gca,'FontWeight','Normal');
    set(gca,'TickLength',[tickLength, tickLength])
    set(gcf,'color','w');
    
    %remove up and right ticks
    %removeURTicks
end