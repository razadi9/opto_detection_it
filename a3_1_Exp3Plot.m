clear
close all

eTitle = true;

load a3_1_Exp3Compute.mat

%%
for iMonkey =1:2


    iX = 80;

    hFig = figure(iMonkey);
    set(hFig,'Units','centimeters')
    if iMonkey ==1 && false
        set(hFig, 'Position',[ 15.4869   13.2997   18.3    6.6])
        hAx2 = subplot(1,2,2);
        hAx1 = subplot(1,2,1);
    else
        set(hFig, 'Position',[ 15.4869   13.2997   18.3/2    6.6])
        hAx2 = [];
        hAx = hAx2;
        clear hAx2
        hAx2 = axes;
    end

    if iMonkey ==1 && false

        hold on
        box off;set(hAx1,'TickDir','out');
        ylabel('Performance','FontSize',9,'FontName','Arial')
        xlabel('Image visiblity'  ,'FontSize',9,'FontName','Arial')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Arial','fontsize',9)
        x = 0:100;
        y = x.^3;
        %y  = y - mean(y)*.7*2;
        y  = y - mean(y)*.7;

        plot(x,y,'LineWidth',1.5)
        a = atand((y(iX+5) -y(iX))./(x(iX+5) -x(iX)) / 20000);
        %a = atand((y(iX+5) -y(iX))./(x(iX+5) -x(iX)) / 20000)*1.2;
        text(x(iX),y(iX)*1.1,{'distortion'},'FontSize',8,'Color',RzTlBx.colorOrder(1) , 'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',a)

        y = -x.^3;
        %y = y-mean(y)*.7*2;
        y = y-mean(y)*.7;
        plot(x,y,'LineWidth',1.5)
        %axis equal
        text(x(iX),y(iX)*1.1,{'hallucination'},'FontSize',8,'Color',RzTlBx.colorOrder(2) , 'HorizontalAlignment','left','VerticalAlignment','top','Rotation',-a)
        yy = y(1);

        y = x*0;
        y  =  yy - y;
        plot(x,y,'LineWidth',1.5)
        %axis equal
        text(x(20),y(20)*1.1,{'hallucination'},'FontSize',8,'Color',RzTlBx.colorOrder(3) , 'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',-a*0)

        if eTitle
            title('Prediction')
            set(get(gca, 'title'),'fontweight','normal')

        end
        y = ylim;
        %ylim(y + diff(y).*[-1 .1])

        xticks(mean(xlim))
        xticklabels({' '})
        yticks(mean(ylim))
        yticklabels({' '})
        %set(gca, 'TickLength',[0 0])

        p = get(gca,'Position');

        axesText = axes('Parent', gcf, ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 1], ...
            'Visible', 'off', ...
            'XLim', [0, 1], ...
            'YLim', [0, 1], ...
            'NextPlot', 'add');

        axes(axesText);
        if eTitle
            text(p(1)-.02,p(2)+p(4)+.01, 'a', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','bold');
        else

            text(p(1)-.05,p(2)+p(4)+.01, 'a.', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','bold');
            text(p(1)-.05+.03,p(2)+p(4)+.01, 'Prediction', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','normal');
        end
    end
    %%

    axes(hAx2)
    hold on





    %set(gca, 'Position', [0.1300+0.7750/2   0.1100    0.7750/2    0.8150])

    nImg = size(img(iMonkey).index,1);
    for iImg = 1:nImg
        ind = img(iMonkey).index(iImg,:);

        y = dPrime(iMonkey).d(ind);
        x = 1:4;

        col = RzTlBx.colorOrder(1);
        col = col+(1-col)*3/4;
        plot(x,y,'-o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',3,'LineWidth',.25)
    end




    nVsb = size(img(iMonkey).index,2);
    y = nan(nVsb+1,1);
    CI= nan(nVsb+1,2);
    for iVsb = 1:nVsb
        ind = img(iMonkey).index(:,iVsb);
        y(iVsb+1) = mean(dPrime(iMonkey).d(ind));
        CI(iVsb+1, 1:2) = prctile( mean(dPrime(iMonkey).Bt(ind,:)), [2.5, 97.5]);
    end
    ind = img(iMonkey).g;
    y(1) = dPrime(iMonkey).d(ind);
    CI(1, 1:2) = prctile( dPrime(iMonkey).Bt(ind,:), [2.5, 97.5]);

    x = 0:4;
    col = RzTlBx.colorOrder(1);
    errorbar(x,y, y - CI(:,1), CI(:,1) - y, '-o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',5,'LineWidth',.5,'CapSize',3)


    xlim([-.2,4.2])
    ylim([min(dPrime(iMonkey).d)-.1,max(dPrime(iMonkey).d)+.1])
    y = ylim;
    ylim([0,y(2)])
    box off;set(hAx2,'TickDir','out');
    xticks(0:4)
    yticks(0:100)
    ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
    xlabel('Image visiblity (a. u.)'  ,'FontSize',9,'FontName','Arial')
    xticklabels({'no image','1','2','3','4'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Arial','fontsize',9)
    xtickangle(0)



    y = [];
    x = [];
    for iVsb = 1:nVsb
        ind = img(iMonkey).index(:,iVsb);
        x = [x;zeros(size(ind))+iVsb];
        y = [y;dPrime(iMonkey).d(ind)'];
    end
    ind = img(iMonkey).g;
    y(end+1) = dPrime(iMonkey).d(ind);
    x(end+1) = 0;



    yl = ylim;
    if eTitle
        if iMonkey == 1
            title('Observation')
            title('Monkey Ph.')
        else
            title('Monkey Sp.')
        end
        set(get(gca, 'title'),'fontweight','normal')

    end

    p = get(gca,'Position');
    xlim([-.5,4.5])
    xl = xlim;
    yl = ylim;
    n = NNN(iMonkey);
    s = num2strComma(n);

    text(xl(2),yl(1)+diff(yl)*.225, sprintf('n = %s trials  ',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')
    if iMonkey ==2
        text(xl(1),yl(1)+diff(yl)*.225, '   Monkey: Sp','FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')
    end

    [pAnova,anovatab,stats] = anova1(y,x,'off');

    [a ,b,c] = multcompare(stats,'Display','off','CType','lsd');
    q = mafdr( a(a(:,1)>=1 & a(:,2)==5,6),'BHFDR', true);
    p = floor(100 * q)/100;
    yl = ylim;
    %ylim('auto')
    for iY = 1:numel(p)
        xx = [iY-.1,5.1]-1;
        yy = [1 1].*max(y)+diff(yl)*.06+ diff(yl).*((4-iY)*.06);
        line(xx, yy, 'Color', 'k')
        text(mean(xx),mean(yy),sprintf('       p=%0.2f       ', p(iY)), 'Color','w','BackgroundColor','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','Helvetica');
        text(mean(xx),mean(yy),sprintf('p=%0.2f', p(iY)), 'Color','k','BackgroundColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',6,'FontName','Helvetica');
    end
    ylim([yl(1), 1.*max(y)+ diff(yl).*(5*.06)])

    [r,p] = corr(x,y,'type','Spearman');
    r = num2str(r,2);
    if p>0.001
        p = ['=',num2str(p,3)];
    else
        p = '< 0.001';
    end
    yl = ylim;
    text(-.3,(yl(1)+yl(2)).*.85, {['r = ' , r],['p ',p]},'FontSize',8,'FontWeight','normal','HorizontalAlignment','left','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    s = {'Monkey Ph', 'Monkey Sp'};
    %text(4,max(ylim),s{iMonkey},'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')


    p = get(gca,'Position');
    if iMonkey ==1
        axesText = axes('Parent', gcf, ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 1], ...
            'Visible', 'off', ...
            'XLim', [0, 1], ...
            'YLim', [0, 1], ...
            'NextPlot', 'add');

        axes(axesText);
        if false
        if eTitle
            text(p(1)-.02,p(2)+p(4)+.01, 'b', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','bold');
        else

            text(p(1)-.05,p(2)+p(4)+.01, 'b.', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','bold');
            text(p(1)-.05+.03,p(2)+p(4)+.01, 'Observation', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial', 'FontWeight','normal');
        end
        end
    end
    for i = 1:5
        a  = axes;
        if iMonkey ==1
            rgbImage = imread(['stim',num2str(i),'.png']);
        else
            rgbImage = imread(['stim2_',num2str(i),'.png']);
        end
        rgbImage = imresize(rgbImage,[220,220]);
        rgbImage = fliplr(rgbImage);
        if i == 5
            rgbImage(:)= 188;
        end
        imshow(rgbImage);
        box off
        axis off
        %%
        if iMonkey ==1 && false
            set(gca, 'Position', [p(3)+.487+(1-i)*.067,p(2)+.025,.1,.1])
        else
            set(gca, 'Position', [p(3)+(1-i)*.154,p(2)+.025,.1,.1])
        end

    end



    %%
    if true
        if iMonkey ==1
            saveas(hFig, 'fig4_Exp3___a3_1.fig')
            saveas(hFig, 'fig4_Exp3___a3_1.pdf')
            saveas(hFig, 'fig4_Exp3___a3_1.png')
        else
            saveas(hFig, 'figED4_Exp3___a3_1.fig')
            saveas(hFig, 'figED4_Exp3___a3_1.pdf')
            saveas(hFig, 'figED4_Exp3___a3_1.png')
        end
    end
end


