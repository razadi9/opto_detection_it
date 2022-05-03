clear
close all

load a1_1_Exp1Compute.mat
e2ndHist = true;


%%
for iMonkey =1


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




    hFig = figure(iMonkey);

    set(hFig,'Units','centimeters')
    set(hFig, 'Position',[ 15.4869   13.2997   18.3   22])


    %% Detection Profile
    hAx1 = axes;
    hold on

    iSite = 2;

    y{40} = [];
    for iImage = 1:40
        dPrime(iMonkey).Bt(iImage,iSite,:);
        c = prctile( dPrime(iMonkey).Bt(iImage,iSite,:),[2.5 97.5]);
        y{iImage} = permute(dPrime(iMonkey).Bt(iImage,iSite,:),[3,1,2]);
        I = y{iImage} <c(1) | y{iImage} >c(2);
        y{iImage}(I) = [];
    end

    %[h,L,MX,MED,bw] = violin(permute(dPrime(iMonkey).Bt(:,iSite,:),[3,1,2]),'edgecolor','k','x',1:40);
    N1 = [1:34,36:40];
    N2 = 35;
    [h,L,MX,MED,bw] = violin(y(N1),'edgecolor','k','x',N1,'facecolor',RzTlBx.colorOrder(1));
    colNoImage = RzTlBx.colorOrder(1)+(1-RzTlBx.colorOrder(1))*.25;
    %colNoImage = [1,1,1]*.5;
    [h,L,MX,MED,bw] = violin(y(N2),'edgecolor','k','x',N2,'facecolor',colNoImage+(1-colNoImage)/2);
    set(h,'EdgeColor',[.5,.5,.5])
    text(N2-.5,dPrime(iMonkey).d(N2,iSite),{'no image'},'Color',colNoImage, 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','top', 'HorizontalAlignment', 'right','FontSize', 8,'FontName', 'helvetica');

    xlim([.5, 40.5])
    ylim([min(cellfun(@min, y)),max(cellfun(@max, y))]+[-.15,.1])
    % h:     figure handle
    % L:     Legend handle
    % MX:    Means of groups
    % MED:   Medians of groups
    % bw:    bandwidth of kernel
    legend off


    plot(1:40, dPrime(iMonkey).d(:,iSite),'k.');

    xticks(10:10:100)
    yticks(0:4:10)
    ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
    xlabel('Image index'  ,'FontSize',9,'FontName','Arial')
    s = {'Monkey Ph', 'Monkey Sp'};
    %text(  20  ,max(ylim),s{iMonkey},'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')





    y = ylim;
    ylim(y + diff(y)/20.*[-1 1])
    title('Detection profile')
    axP=[0.05    0.05   0.9   .95 ];
    h = axP(4)*1/3-.075;
    w = axP(3)*4/5-.025;
    axesMain( [axP(1)    1-(axP(2)+h)   w   h ])
    x= xlim;
    y= xlim;

    n = NNN(iMonkey).(['stim',num2str(iSite)])+NNN(iMonkey).noStim;
    s = num2strComma(n);

    text(min(xlim)+diff(xlim)*.03*0,min(ylim)+diff(ylim)*.05*0, sprintf('  n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

    p = get(gca,'Position');

    axesText = axes('Parent', gcf, ...
        'Units', 'normalized', ...
        'Position', [0, 0, 1, 1], ...
        'Visible', 'off', ...
        'XLim', [0, 1], ...
        'YLim', [0, 1], ...
        'NextPlot', 'add');

    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'a', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');


    %% Detection Profile - hist

    hAx2 = axes;
    hold on

    dPerm = permute(std(dPrime(iMonkey).dPerm(:,iSite,:),0,1),[3,2,1]);
    dObs = std(dPrime(iMonkey).d(:,iSite,:),0,1);


    histPerm(dPerm, dObs, true, true)
    yticks(.10:.7:10)
    ylim([.1,.8])


    xh = get(gca,'ylabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(1) = p(1)*1.5 ;        % double the distance, 
                       % negative values put the label below the axis
    set(xh,'position',p);   % set the new position


    ylabel('Performance SD (\sigma_{d''})','FontSize',9,'FontName','Arial')
    title('Permutation Test','FontSize',9,'FontName','Arial')


    p = get(hAx1,'Position');
    h = p(4);
    w = axP(3)*1/5-.025;

    axesMain( [axP(1)+axP(3)*4/5    p(2)   w   h ])


    
    p = get(gca,'Position');
    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'b', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');


    
    %% Correlation Between and Within
    clear dPrime
    load('a1_2_Exp1ComputeCorr.mat')
    hAx3 = axes;
    hold on



    x1 = [dPrime(iMonkey).rWithin1; dPrime(iMonkey).rWithin2];
    x2 = [dPrime(iMonkey).rBetween1;dPrime(iMonkey).rBetween2];


    [f,x]= ecdf(x1);
    line(x,f,'Color',RzTlBx.colorOrder(4))
    [f,x]= ecdf(x2);
    line(x,f,'Color',RzTlBx.colorOrder(5)-[.2,0,.1])


    %xticks(0:.2:1);
    yticks(0:.5:1)

    %plot([1,5],(yl(1)-.05)*[1,1],'k-')
    ylabel({'Cumulative probability'},'FontSize',9,'FontName','Arial')
    yticks([0,1])
    xlabel('Correlation coefficient','FontSize',9,'FontName','Arial')
    title('Correlation of stimulation sites')

    xlim(prctile([x1;x2],[0,100]))
    %xlim([.3,.9])
    xticks([])
    xticks(0:.2:1)


    [y,x] = ecdf(x1);
    [~,I] = min(abs(y-.55));
    [~,J] = min(abs(y-.92));
    a = atand((y(J)-y(I)+.05)./(x(J)-x(I))/(diff(ylim)/diff(xlim)));
    text(x(I)+.01,y(I),'Within sites','FontSize',8, 'Color',RzTlBx.colorOrder(4),'HorizontalAlignment','left','VerticalAlignment','top', 'Rotation',a);

    [y,x] = ecdf(x2);
    [~,I] = min(abs(y-.55));
    [~,J] = min(abs(y-.92));
    a = atand(((y(J)-y(I))./(x(J)-x(I))/(diff(ylim)/diff(xlim))));
    text(x(I)-.01,y(I),'Between sites','FontSize',8, 'Color',RzTlBx.colorOrder(5)-[.2 0 .1],'HorizontalAlignment','left','VerticalAlignment','bottom', 'Rotation',a);


    [y1,x1] = ecdf(x1);
    [y2,x2] = ecdf(x2);
    [~,I1] = min(abs(y1-.4));
    [~,I2] = min(abs(y2-.4));
    %plot([x1(I1),x2(I2)]+[-.01 0.01],[y1(I1),y2(I2)],'k');
    xx1 = x1(I1);
    xx2 = x2(I2);
    text(mean([xx1, xx2]),.4,'p < 0.001' ,'FontSize',8,'Color','k','BackgroundColor','none','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',0);


    [~,I1] = min(abs(y1-.5));
    [~,I2] = min(abs(y2-.5));
    xx1 = x1(I1);
    xx2 = x2(I2);



    %text(mean([x2(I2)]),mean([y1(I1),y2(I2)]),['\leftarrow'] ,'FontSize',18,'Color','k','BackgroundColor','none','HorizontalAlignment','left','VerticalAlignment','middle','Rotation',0,'FontWeight','bold');
    %text(mean([x1(I1)]),mean([y1(I1),y2(I2)]),['\rightarrow'] ,'FontSize',18,'Color','k','BackgroundColor','none','HorizontalAlignment','right','VerticalAlignment','middle','Rotation',0,'FontWeight','bold');
    %stext(mean([x1(I1)]),mean([y1(I1),y2(I2)]),['\lefttarrow'] ,'FontSize',15,'Color','k','BackgroundColor','none','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',0,'FontWeight','bold');
    %plot([1,1]*x1(I1),[1,0]*y1(I1),'Color',RzTlBx.colorOrder(1))
    %plot([1,1]*x2(I2),[1,0]*y2(I2),'Color',RzTlBx.colorOrder(2))

%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))

    xticks(0:.4:1)

    h = axP(4)*.8/3-.05;
    w = axP(3)*2/5-.025;
    p = get(hAx1,'Position');
    axesMain( [axP(1)    p(3)-.1-h  w   h ])


    [xx yy] = datc2figc([xx1, xx2],[.5, .5]);
    ar = annotation('doublearrow',xx,yy,'Head1Length',7,'Head1Width',7,'Head2Length',7,'Head2Width',7);
    
    p = get(gca,'Position');
    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'c', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');
    




        a  = axes;
        rgbImage = imread('array_schematic_beveled copy.png','BackgroundColor',[1 1 1]);
        rgbImage = imresize(rgbImage,.5);
        imshow(rgbImage);
        set(gca,'position',[p(1)+p(3)*.02,p(2)+p(4)*.25,.125,.125])
        box off
        axis off
    
    1==1
    %% Correlation Between and Within - hist
    %{
    hAx4 = axes;
    hold on

    dPerm = dPrime(iMonkey).rPrm(2:end);
    dObs = dPrime(iMonkey).rPrm(1);

    
    histPerm(dPerm, dObs, false, false)
    ylim([.9 1])
    ylabel('Correlation coefficient','FontSize',7,'FontName','Arial')
    title('Permutation Test','FontSize',7,'FontName','Arial')
    yticks(0:.1:1)

    p = get(hAx3,'Position');
    h = p(4);
    w = axP(3)*2/5-.025;

    axesMain( [axP(1)+axP(3)*3/5    p(2)   w   h ])


    %}
    %% Psychometric
    iSite   = 2;
    load('a2_1_Exp2ComputeInt.mat')
    hAx4 = axes;
    p = get(hAx3,'Position');
    h = p(4);
    if e2ndHist
        w = axP(3)*2/5-.025*2;
    else
        w = axP(3)*3/5-.025*2;
    end
    hold on
    xlim([-5,105])

    axesMain( [axP(1)+axP(3)*2/5+.025    p(2)   w   h ])

    if iMonkey ==1
        N1 = 1:4;
        N2 = 5;
    else
        N1 = [1,2,3,5];
        N2 = 4;
    end
    h = plot(UIntAll{iMonkey}*1000, squeeze(dPrime(iMonkey).d(N1,iSite,:) )','-','Marker','o','MarkerSize',5,'LineWidth',1);
    for iH= 1:numel(h)
        set(h(iH),'color',RzTlBx.colorOrder(iH+1))
    end
    plot(UIntAll{iMonkey}*1000, squeeze(dPrime(iMonkey).d(N2,iSite,:) )','-','Marker','s','MarkerSize',5,'LineWidth',1,'color',colNoImage)
    x = UIntAll{iMonkey}(4:5)*1000;
    y= squeeze(dPrime(iMonkey).d(N2,iSite,4:5) );

    a  = atand((diff(y)./diff(x))./(diff(ylim)./diff(xlim)))-10;
    text(mean(x), mean(y)-.2,{'no image'},'Color',colNoImage, 'BackgroundColor', 'none','Rotation', a,'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 7,'FontName', 'helvetica');
    xlabel('Illumination power (\muW)','FontSize',9,'FontName','Arial')
    ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
    %title('Permutation Test','FontSize',7,'FontName','Arial')
    yticks(-2:2:10)
    %xticks(0:50:1000)
    title('Intensity psychometric functions')
    n = NNN(iMonkey).noStim + NNN(iMonkey).(['stim',num2str(iSite)]);
    y = ylim;
    ylim([-.5,4])
    s = num2strComma(n);
    
    text(max(xlim)-diff(xlim)*.05*0,min(ylim)+diff(ylim)*.00, sprintf('n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')


    1==1
    %h = axP(4)*1.2/3-.05;
    %w = axP(3)*3.5/5-.025;
    %p = get(hAx3,'Position');
    %axesMain( [axP(1)    axP(1)  w   h ])

    p = get(gca,'Position');
    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'd', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');

    %%
    %% Correlation Between and Within - hist
    %
    if e2ndHist
    hAx6 = axes;
    hold on

    iSite = 2;
    dPerm = dPrime(iMonkey).rPerm(iSite,:);
    dObs = dPrime(iMonkey).rPerm(iSite,1);

    histPerm(dPerm, dObs, true, false)

    ylabel('Corr. coefficient mean','FontSize',7,'FontName','Arial')
        xh = get(gca,'ylabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(1) = p(1)*1.5 ;        % double the distance, 
                       % negative values put the label below the axis
    set(xh,'position',p);   % set the new position

    title('Permutation Test','FontSize',7,'FontName','Arial')

    y = ylim;
    ylim(y+diff(y)*.1.*[0 1])
    ylim([-1 1])
    yticks([-1,1])
    


    
    

    p = get(hAx4,'Position');
    p1 = p(1)+p(3)+.025;

    axesMain( [p1    p(2)   1-p1-.025*3  p(4) ])
    end

    %}

    saveas(hFig, 'fig3_Exp1_2___a1_9.fig')
    saveas(hFig, 'fig3_Exp1_2___a1_9.pdf')
    saveas(hFig, 'fig3_Exp1_2___a1_9.png')
    fakdjfa
end




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
    text(mean(xl),  prctile(dPerm, 98),sprintf('Null distribution'),'Color',RzTlBx.colorOrder(2), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','bottom', 'HorizontalAlignment', 'center','FontSize', 8,'FontName', 'helvetica')

    if dObsTextDwon
        text(mean(xl), dObs-diff(yl)*.05 ,sprintf('Observed'),'Color',RzTlBx.colorOrder(1), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','cap', 'HorizontalAlignment', 'center','FontSize', 8,'FontName', 'helvetica')
    else
        text(mean(xl), dObs+diff(yl)*.05 ,sprintf('Observed'),'Color',RzTlBx.colorOrder(1), 'BackgroundColor', 'none','Rotation', 0,'VerticalAlignment','baseline', 'HorizontalAlignment', 'center','FontSize', 8,'FontName', 'helvetica')
    end
end
xlabel('Frequency','FontSize',9,'FontName','Arial')
xl = xl+diff(xl)*[0,.05];
xlim(xl)
a = get(gca,'XTickLabel');

end


