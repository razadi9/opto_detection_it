clear
close all

load a0_1_trainingPhaseCompute2
cpSz = 0;
%%
for iMonkey =1

    hFig = figure(iMonkey);

    set(hFig,'Units','centimeters')
    set(hFig, 'Position',[ 15.4869   13.2997   8.8   6.6])
    hAx1 = axes; 
    hold on



    x = 1:size(stimReport(iMonkey).m,2);

    xDrift =[-.15,0,.15];
    for iTrialType = 3:-1:1
        y = stimReport(iMonkey).m(iTrialType,:);
        
        CI1 = stimReport(iMonkey).btCI(iTrialType,:,1);
        CI2 = stimReport(iMonkey).btCI(iTrialType,:,2);

        if iMonkey ==1 && iTrialType == 3
            y(1:11) = nan;
            CI1(1:11) = nan;
            CI2(1:11) = nan;
        end

        col = RzTlBx.colorOrder(iTrialType);

        y(end) = nan; % this is the average of three last sessions
        I = ~isnan(y);
        
        if false
            line(x(I), y(I),'color',col)
            fill([x(I),flip(x(I))], [CI2(I),flip(CI1(I))], col, 'EdgeColor','none','FaceAlpha',.5)
        else
            errorbar(x(I)+xDrift(iTrialType), y(I), y(I)-CI1(I) , CI1(I)-y(I) ,'color',col, 'CapSize',cpSz,'LineWidth',.5)
        end

        if iTrialType  ==2
            %x1stSigSession = find(permTest(iMonkey).p(:,1)<.05,1,'first');
            p = [stimReport(iMonkey).chiSq(:,1).p]';
            x1stSigSession = 4;
            text(x1stSigSession,CI2(x1stSigSession)+.07,{sprintf('     p = %0.2f ',p(x1stSigSession,1)),' '},'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','helvetica')
            text(x1stSigSession,CI2(x1stSigSession)+.05,{['   '],'\downarrow'},'FontSize',9,'FontWeight','bold'    ,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
    end


    
    
    
    plot([1,1]*39.5,[.05,.95],'k:','LineWidth',.33)
    plot([1,1]*42.5,[.05,.95],'k:','LineWidth',.33)

    xlim([0,x(end)+1])
    ylim([0,1])

    box off;set(hAx1,'TickDir','out');
    xticks(10:10:1000)
    yticks(0:1)
    ylabel('Stimulation report rate','FontSize',9,'FontName','Arial')
    xlabel('Session number','FontSize',9,'FontName','Arial')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Arial','fontsize',9)
    %title(sprintf('Training phase'))
  


    p = [0.1300    0.1440    0.7750   0.7810]; % defult
    set(hAx1,'Position', p + [0 , 0 , -.1 , 0])
    p = get(hAx1,'Position');
    set(hAx1, 'Position', p .*[1,1,1,.95]+[0,p(4).*(1-.95),0,0])

    
    
    n = sum(cellfun(@numel,N_trialTypeDay{iMonkey}( 1:end-1,:)),'all');
    s = num2strComma(n);
    text(1,0, sprintf('n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')


    text(20,.9,'stimulation','FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none', 'Color',RzTlBx.colorOrder(2))
    text(1,.15,'non-stimulation','FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none', 'Color',RzTlBx.colorOrder(1))
    text(20,.03,'catch','FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none', 'Color',RzTlBx.colorOrder(3))

    
    % violin plot
    hAx2 = axes;
    hold on
    p = get(hAx1,'Position');
    set(hAx2,'Position', p + [p(3),0,-p(3)+.1+.05,0])
    p = get(gcf,'Position');
    set(gcf, 'Position', p .*[1,1,1,.8])
    

    a = get(gca, 'XAxis');
    set(a, 'Visible','off');
    a = get(gca, 'YAxis');
    set(a, 'Visible','off');
    hold on
    x = [1,1,1.7];

    
    for iTrialType = 1:3
        y = stimReport(iMonkey).m(iTrialType,end);
        CI1 = stimReport(iMonkey).btCI(iTrialType,end,1);
        CI2 = stimReport(iMonkey).btCI(iTrialType,end,2);
        col = RzTlBx.colorOrder(iTrialType);
        if false
            errorbar(x(iTrialType),y,y-CI1,CI2-y,'Color', col,'Marker','o','MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',2,'CapSize',cpSz,'LineWidth',1)
        elseif false
            fill([-1 -1 1 1]/2.5+x(iTrialType),[CI1,CI2,CI2,CI1], col,'FaceAlpha',.5,'EdgeColor',col,'LineWidth',.25,'LineStyle','none')
            line([-1 1]/3+x(iTrialType),[y,y], 'color',    col, 'linewidth',1)
        else
            y = squeeze(stimReport(iMonkey).btData(iTrialType,end,:));
            yci = prctile(y,[2.5,97.5]);
            y = y(y>=yci(1) & y<=yci(2));
            violin(y,'x', x(iTrialType),'facecolor', col)
            plot( x(iTrialType),stimReport(iMonkey).m(iTrialType,end),'k.')

        end
    end
    xlim([.5,2])
    ylim([0,1])

    xx = [-1 -1]/2.5+x(1) -.1;
    yy = stimReport(iMonkey).m([1,2],end);
    line(xx,yy, 'color', 'k', 'linewidth',.5)
    line(xx(1)+[0,.1],yy(1).*[1,1], 'color', 'k', 'linewidth',.5)
    line(xx(1)+[0,.1],yy(2).*[1,1], 'color', 'k', 'linewidth',.5)
    line(xx(1)+[0,.1],mean(yy).*[1,1], 'color', 'k', 'linewidth',.5)
    text(xx(1)+.1,mean(yy),' p < 0.001','FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    
    
    xx = x([1,3]);
    yy = [1,1]*min(y)-.03;
    line(xx,yy-.02, 'color', 'k', 'linewidth',.5)
    line(xx(1).*[1,1],yy+[0,-.02], 'color', 'k', 'linewidth',.5)
    line(xx(end).*[1,1],yy+[0,-.02], 'color', 'k', 'linewidth',.5)
    line(mean(xx).*[1,1].*[1,1],yy+[-.02,-.035], 'color', 'k', 'linewidth',.5)
         

    p  = stimReport(iMonkey).chiSq(end,1).p;
    
    
    text(mean(xx),mean(yy),['n.s.'],'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','top','FontName','helvetica','Rotation',0, 'BackgroundColor','none')
    
    %s = {'Monkey Ph', 'Monkey Sp'};
    %text(.4,1.04,s{iMonkey},'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    
if iMonkey ==1
     [stimReport(iMonkey).chiSq([4,end],1).p]
     [stimReport(iMonkey).chiSq([4,end],1).chi2]
     [sum(stimReport(iMonkey).chiSq(4,1).tbl(:)), sum(stimReport(iMonkey).chiSq(end,1).tbl(:))]
     ch2df = (2-1)*(2-1)

     % all compariosns for catch and stimulation trials.
     min([stimReport(iMonkey).chiSq([17:end-1],2).p])
elseif iMonkey == 2
    find(([stimReport(iMonkey).chiSq(19:end,1).p])'<.05)
    sum([stimReport(iMonkey).chiSq(19+11,1).tbl]','all')
    ([stimReport(iMonkey).chiSq(19+11,1).p])
    min([stimReport(iMonkey).chiSq(19:end,2).p])
    
end


    saveas(hFig, 'fig2_TrainingPhase___a0_2.fig')
    saveas(hFig, 'fig2_TrainingPhase___a0_2.pdf')
    saveas(hFig, 'fig2_TrainingPhase___a0_2.png')
lkajdfalj
end


