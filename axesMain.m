
function axesMain(p)
    set(gca,'Position',p)
    hold on
    set(gca,'TickDirMode','manual')
    box off;set(gca,'TickDir','out');
    set(gca,'FontName','Arial','fontsize',9)
end
