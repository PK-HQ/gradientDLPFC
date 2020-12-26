function removeURTicks()
%REMOVE UP AND RIGHT AXES AND TICKS
set(gca,'box','off')

isholdonque = ishold;
hold on
ax = axis;
%plot(ax(2)*[1,1],ax(3:4),'k','linewidth',2.5)
%plot(ax(1:2),ax(4)*[1,1],'k','linewidth',2.5)
if isholdonque == 0
    hold off
end
end