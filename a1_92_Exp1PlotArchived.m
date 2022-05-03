clear
close all

load a1_1_Exp1Compute.mat



%%
for iMonkey =1:2


    for iSite = 1:size(dPrime(1).dPerm,2)
        for iImage = 1:40
            I = isnan(dPrime(iMonkey).dPerm(iImage,iSite,:));
            if any(I)
                warning('there is NaN value in the data.')
                N = find(~I);
                N = N(randperm(numel(N)));
                N = datasample(N,numel(I),'Replace',true);
                dPrime(iMonkey).dPerm(iImage,iSite,:) = N;
            end
        end
    end

    % this is because sp data has an NaN on one of the control LEDs
    I = isnan(dPrime(iMonkey).d);
    dPrime(iMonkey).d(I) = 0;

    d = dPrime(iMonkey).d;
    dp = dPrime(iMonkey).dPerm;

    pValue{iMonkey} = sum(std(d,1)  < squeeze(std(dp,1,1))')./size(dPrime(iMonkey).dPerm,3)

    d = [mean(dPrime(iMonkey).d([1:34,36:end],:,:));(dPrime(iMonkey).d([35],:,:))] ;
    dp = [mean(dPrime(iMonkey).dPerm([1:34,36:end],:,:));(dPrime(iMonkey).dPerm([35],:,:))] ;


    % analyze of gray
    pValueGray{iMonkey} = sum(std(d,1)  < squeeze(std(dp,1,1))')./size(dPrime(iMonkey).dPerm,3)
    pValueGray{iMonkey}  = mafdr(pValueGray{iMonkey},'BH',true)

end


hFig = figure;

set(hFig,'Units','centimeters')
set(hFig, 'Position',[ 15.4869   13.2997   18.3   22])
axesText = axes('Parent', gcf, ...
    'Units', 'normalized', ...
    'Position', [0, 0, 1, 1], ...
    'Visible', 'off', ...
    'XLim', [0, 1], ...
    'YLim', [0, 1], ...
    'NextPlot', 'add');


%% Detection Profile
for iMonkey = 1:2
    if iMonkey ==1
        hAx1 = axes;
        hold on

        axP=[0.05    0.05   0.9   .95 ];
        h = axP(4)*1/3-.075;
        w = axP(3)*4/5-.025;
        w = axP(3);
        axesMain( [axP(1)    1-(axP(2)+h)   w   h ])
    else
        hAx2 = axes;
        hold on
        h = axP(4)*1/3-.075 ;
        w = axP(3);
        p = get(hAx1,'Position');
        axesMain( [axP(1)    p(2)-.1-h  w   h ])



    end
    if iMonkey ==2
        x = dPrime(iMonkey).Bt   (:,1,:);dPrime(iMonkey).Bt   (:,1,:) = dPrime(iMonkey).Bt   (:,3,:);dPrime(iMonkey).Bt   (:,3,:) = x;
        x = dPrime(iMonkey).CI   (:,1,:);dPrime(iMonkey).CI   (:,1,:) = dPrime(iMonkey).CI   (:,3,:);dPrime(iMonkey).CI   (:,3,:) = x;
        x = dPrime(iMonkey).d    (:,1  );dPrime(iMonkey).d    (:,1  ) = dPrime(iMonkey).d    (:,3  );dPrime(iMonkey).d    (:,3  ) = x;
        x = dPrime(iMonkey).dPerm(:,1,:);dPrime(iMonkey).dPerm(:,1,:) = dPrime(iMonkey).dPerm(:,3,:);dPrime(iMonkey).dPerm(:,3,:) = x;


    end
    for iSite = 1:2

        clear y
        y{40} = [];
        for iImage = 1:40
            dPrime(iMonkey).Bt(iImage,iSite,:);
            c = prctile( dPrime(iMonkey).Bt(iImage,iSite,:),[2.5 97.5]);
            y{iImage} = permute(dPrime(iMonkey).Bt(iImage,iSite,:),[3,1,2]);
            I = y{iImage} <c(1) | y{iImage} >c(2);
            y{iImage}(I) = [];
        end

        N1 = [1:34,36:40];
        N2 = 35;
        k = 2;
        x = (k:k:40*k)+ (iSite-.5)*.75 ;
        [h,L,MX,MED,bw] = violin(y(N1),'edgecolor','k','x',x(N1),'facecolor',RzTlBx.colorOrder(iMonkey));
        set(h,'EdgeColor','none')
        [h,L,MX,MED,bw] = violin(y(N2),'edgecolor','k','x',x(N2),'facecolor',RzTlBx.colorOrder(iMonkey));
        set(h,'EdgeColor',RzTlBx.colorOrder(iMonkey))
        set(h,'FaceAlpha',1)
        plot(x, dPrime(iMonkey).d(:,iSite),'.','Color',RzTlBx.colorOrder(iMonkey));
        set(h,'FaceAlpha',.1)
        set(h,'EdgeColor',RzTlBx.colorOrder(iMonkey))
        %text(N2+.5,dPrime(iMonkey).d(N2,iSite),{'gray background'},'Color',(RzTlBx.colorOrder(1)+.2)/1.2, 'BackgroundColor', 'none','Rotation', 90,'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 7,'FontName', 'helvetica');

        %xlim([.5, 40.5])
        % h:     figure handle
        % L:     Legend handle
        % MX:    Means of groups
        % MED:   Medians of groups
        % bw:    bandwidth of kernel
        legend off
    end
    xlim([x(1),x(end)]+[-k 1])

    ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
    xlabel('Image number'  ,'FontSize',9,'FontName','Arial')
    %s = {'Monkey Ph', 'Monkey Sp'};
    %text(  20  ,max(ylim),s{iMonkey},'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')


    y = ([min(min(min(dPrime(iMonkey).CI(:,1:2,:)))), max(max(max(dPrime(iMonkey).CI(:,1:2,:)))) ]);
    y = y+[-1 1].*diff(y)/10;

    ylim(y);
    


    title('Detection profile')

    %n = NNN(iMonkey).(['stim',num2str(iSite)])+NNN(iMonkey).noStim;
    %text(min(xlim)+diff(xlim)*.03,min(ylim)+diff(ylim)*.05, sprintf('n = %d trials',n),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    p = get(gca,'Position');

    x1 = (k:k:40*k)+ (1-.5)*.75 ;
    x2 = (k:k:40*k)+ (2-.5)*.75 ;
    plot([x1',x2']', dPrime(iMonkey).d(:,1:2)','-','Color','k');

    x = xlim;
    xticks(prctile( x(1):x(end), [25, 50, 75]))
    xticklabels({' '})
    yticks(-2:2:10)

    s = {'Monkey Ph', 'Monkey Sp'};
    text(  mean(xlim)  ,max(ylim),s{iMonkey},'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','top','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    if iMonkey ==1

        axP=[0.05    0.05   0.9   .95 ];
        h = axP(4)*1/3-.075;
        w = axP(3)*4/5-.025;
        w = axP(3);
        axesMain( [axP(1)    1-(axP(2)+h)   w   h ])

        p = get(gca,'Position');
        
        axes(axesText);

        text(p(1)-.02,p(2)+p(4)+.02, 'a', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial');
        
    else
        h = axP(4)*1/3-.075 ;
        w = axP(3);
        p = get(hAx1,'Position');
        axesMain( [axP(1)    p(2)-.1-h  w   h ])


        p = get(gca,'Position');
        axes(axesText);
        text(p(1)-.02,p(2)+p(4)+.02, 'b', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial');

    end


    
end


saveas(hFig, 'figExt2_Exp1_2___a1_92.fig')
saveas(hFig, 'figExt2_Exp1_2___a1_92.pdf')
saveas(hFig, 'figExt2_Exp1_2___a1_92.png')





function histPerm(dPerm, dObs, dObsTextDwon, eText)

h = histogram(dPerm,'orientation','horizontal','EdgeColor','none','FaceColor',RzTlBx.colorOrder(2));
set(get(h,'Parent'),'xdir','r')
set(gca, 'YAxisLocation','right')

xticks([]);
yticks(0:.3:10)
xl = xlim;
yl = ylim;
plot(xl,dObs*[1,1],'color',RzTlBx.colorOrder(1),'LineWidth',1)

if eText
    text(mean(xl),  max(dPerm),sprintf('Null distribution\nRandomly labeled (n = 1e4)'),'Color',RzTlBx.colorOrder(2), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','bottom', 'HorizontalAlignment', 'center','FontSize', 6,'FontName', 'helvetica')

    if dObsTextDwon
        text(mean(xl), dObs-diff(yl)*.05 ,sprintf('Observed'),'Color',RzTlBx.colorOrder(1), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','cap', 'HorizontalAlignment', 'center','FontSize', 6,'FontName', 'helvetica')
    else
        text(mean(xl), dObs+diff(yl)*.05 ,sprintf('Observed'),'Color',RzTlBx.colorOrder(1), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','baseline', 'HorizontalAlignment', 'center','FontSize', 6,'FontName', 'helvetica')
    end
end
xlabel('Frequency','FontSize',9,'FontName','Arial')
xl = xl+diff(xl)*[0,.05];
xlim(xl)
a = get(gca,'XTickLabel');

end
