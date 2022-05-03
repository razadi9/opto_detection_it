clear
close all

load a1_1_Exp1Compute.mat


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
    moveyLabel(.1);
    set(get(gca, 'title'),'fontweight','normal')
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
    y = ylim;
    ylim([0,y(2)])
    yticks([0,fix(y(2)*10)/10])


    moveyLabel(1.5);


    ylabel('d'' SD','FontSize',9,'FontName','Arial')
    %title('Permutation Test','FontSize',9,'FontName','Arial')


    p = get(hAx1,'Position');
    h = p(4);
    w = axP(3)*1/5-.025;

    axesMain( [axP(1)+axP(3)*4/5    p(2)   w   h ])
    set(gca,'TickLength',get(gca,'TickLength')*3)






    %% Correlation Between and Within
    clear dPrime
    load('a1_2_Exp1ComputeCorr.mat')
    dPrime(iMonkey).rWithin1 = percent95(dPrime(iMonkey).rWithin1);    
    dPrime(iMonkey).rWithin2 = percent95(dPrime(iMonkey).rWithin2);    
    dPrime(iMonkey).rBetween = percent95(dPrime(iMonkey).rBetween);    
    
    %%
    hAx3 = axes;
    hold on
    fun = @(x)   x(1:fix(numel(x)/1e4):fix(size(x)));
    if numel(dPrime(iMonkey).rWithin1)>1e4
        [h,L,MX,MED,bw, medianLine] = violin( {fun(dPrime(iMonkey).rWithin1), fun(dPrime(iMonkey).rWithin2)}, 'plotmedian', true, 'x',[1,2],'facecolor',RzTlBx.colorOrder(2));
        set(medianLine,'color',[1,1,1]*.5)

        [h,L,MX,MED,bw, medianLine] = violin( {fun(dPrime(iMonkey).rBetween)}                               , 'plotmedian', true, 'x',[4], 'facecolor', RzTlBx.colorOrder(1));
        set(medianLine,'color',[1,1,1]*.5)

    else
        [h,L,MX,MED,bw, medianLine] = violin( {dPrime(iMonkey).rWithin1, dPrime(iMonkey).rWithin2, dPrime(iMonkey).rBetween}, 'plotmedian', true, 'x',[1,2,4]);
    end

    legend off
    xticks([1,2,4]);
    row1 = {'Within' 'Within' 'Between'};
    row2 = {' Site1' ' Site2' '   Sites' };
    labelArray = [row1; row2; ];
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    hAx3.XTickLabel = tickLabels;

    xtickangle(0)

    xlim([.5 , 4.5]);
    %y = ylim;
    %ylim(y + diff(y)/20.*[-1 1])
    ylim([fix((min(dPrime(iMonkey).rBetween))*10)/10-.1,1])
    yticks(ylim)
    yyy = ylim;

    %plot([1,5],(yl(1)-.05)*[1,1],'k-')
    ylabel('Correlation coefficient','FontSize',7,'FontName','Arial')
    %title({'Correlation of', 'stimulation sites'})
    
    h = axP(4)*.8/3-.05;
    w = axP(3)*1.35/5-.025;
    p = get(hAx1,'Position');
    axesMain( [axP(1)    p(3)-.1-h  w   h ])


    s = num2strComma(nAllTrials(iMonkey));

        text(min(xlim)+diff(xlim)*.03*0,min(ylim)+diff(ylim)*.05*0, sprintf('  n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

        moveyLabel(-1.5);
    set(gca,'TickLength',get(gca,'TickLength')*3)
    p = get(gca,'Position');
    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'b', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');
    

        p = get(hAx3, 'Position');
        a  = axes;
        rgbImage = imread('array_schematic_beveled copy.png','BackgroundColor',[1 1 1]);
        rgbImage = fliplr(rgbImage);
        rgbImage = imresize(rgbImage,.5);
        %rgbImage = fliplr(rgbImage);
        imshow(rgbImage);
        set(gca,'position',[p(1)+p(3)*.05, p(2)+p(4)*(-.05), .1, .1])
        box off
        axis off


1==1

    %% Correlation Between and Within - hist
    hAx4 = axes;
    hold on

    dPerm = dPrime(iMonkey).rPrm(2:end);
    dObs = dPrime(iMonkey).rPrm(1);


    histPerm(dPerm, dObs, false, false)
    ylim(yyy)
    yticks(yyy)
    ylabel('Correlation coefficient','FontSize',7,'FontName','Arial')
    %title('Permutation Test','FontSize',7,'FontName','Arial')
    set(gca,'TickLength',get(gca,'TickLength')*3)
    moveyLabel(1.5);


    p = get(hAx3,'Position');
    h = p(4);
    w = axP(3)*.8/5-.025;

    axesMain( [p(1)+p(3)+.025, p(2),   w   h ])





    %% Psychometric
    iSite   = 2;
    load('a2_1_Exp2ComputeIntFit.mat')
    
    min(reshape(dPrime(2).r2(:,1,:),[],1))
    max(reshape(dPrime(2).r2(:,1,:),[],1))
    mean(reshape(dPrime(2).r2(:,1,:),[],1))


    %iMonkey =1;sum(dPrime(iMonkey).s(:,1) < dPrime(iMonkey).s,2)/size(dPrime(iMonkey).s,2) 
    %iMonkey =2;sum(dPrime(iMonkey).s(:,1) < dPrime(iMonkey).s,2)/size(dPrime(iMonkey).s,2) 

    hAx5 = axes;
    p = get(hAx4,'Position');
    h = p(4);
    w = axP(3)*1.65/5-.025;
    hold on
    xlim([-5,105])


    if iMonkey ==1
        N1 = 1:4;
        N2 = 5;
    else
        N1 = [1,2,3,5];
        N2 = 4;
    end

    for iImage= 1:5
        if iImage == N2
            col = colNoImage;
            m = 's';

        else
            col = RzTlBx.colorOrder(iImage+1);
            m = 'o';
        end
        y = squeeze(dPrime(iMonkey).d(iSite,iImage,:) )';

        eFitPolot = false;
        if eFitPolot
            plot(UIntAll{iMonkey}*1000, y,'--','Marker',m,'MarkerSize',5,'LineWidth',.5, 'Color',col);
            x = dPrime(iMonkey).fita( iSite, 1    )*UIntAll{iMonkey}.^dPrime(iMonkey).fitb( iSite, 1    );
            y = x* dPrime(iMonkey).a(iSite, 1,iImage);
            plot(UIntAll{iMonkey}*1000, y,'-','Marker','none','MarkerSize',5,'LineWidth',1, 'Color',col);
            if iImage == N2
                x = UIntAll{iMonkey}(4:5)*1000;
                y= squeeze(dPrime(iMonkey).d(iSite,N2,4:5) );
                a  = atand((diff(y)./diff(x))./(diff(ylim)./diff(xlim)))-10;
                text(mean(x), mean(y)-.2,{'no image'},'Color',colNoImage, 'BackgroundColor', 'none','Rotation', a,'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 7,'FontName', 'helvetica');
            end
        
        else
            plot(UIntAll{iMonkey}*1000, y,'-','Marker',m,'MarkerSize',5,'LineWidth',1, 'Color',col);
            if iImage == N2
                x = UIntAll{iMonkey}(4:5)*1000;
                y= squeeze(dPrime(iMonkey).d(iSite,N2,4:5) );
                a  = atand((diff(y)./diff(x))./(diff(ylim)./diff(xlim)))-10;
                text(mean(x), mean(y)-.2,{'no image'},'Color',colNoImage, 'BackgroundColor', 'none','Rotation', a,'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 7,'FontName', 'helvetica');
            end
        end

    end

    xlabel('Illumination power (\muW)','FontSize',9,'FontName','Arial')
    ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
    %title('Permutation Test','FontSize',7,'FontName','Arial')
    yticks(-4:4:10)
    %xticks(0:50:1000)
    %title({'Intensity', 'psychometric functions'})
    n = NNN(iMonkey).noStim + NNN(iMonkey).(['stim',num2str(iSite)]);
    y = ylim;
    ylim([-.5,4])
    s = num2strComma(n);
    axesMain([p(1)+p(3)+.025*4, p(2), w, h])

    text(max(xlim)-diff(xlim)*.05*0,min(ylim)+diff(ylim)*.00, sprintf('n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')


    1==1
    %h = axP(4)*1.2/3-.05;
    %w = axP(3)*3.5/5-.025;
    %p = get(hAx3,'Position');
    %axesMain( [axP(1)    axP(1)  w   h ])
    set(gca,'TickLength',get(gca,'TickLength')*3)    
    moveyLabel(.75);

    p = get(gca,'Position');
    axes(axesText);
    text(p(1)-.02,p(2)+p(4)+.02, 'c', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');

    %% Psychometric - hist
    hAx6 = axes;
    hold on

    histPerm(dPrime(iMonkey).s(iSite,2:end), dPrime(iMonkey).s(iSite,1), false, false)
    y = ylim.*[0,1];
    ylim(y)
    yticks(y)
    ylabel('Coefficient SD','FontSize',7,'FontName','Arial')
    set(gca,'TickLength',get(gca,'TickLength')*3)
    
    moveyLabel(1.5);


    p = get(hAx5,'Position');
    h = p(4);
    w = axP(3)*.8/5-.025;

    axesMain( [p(1)+p(3)+.025, p(2),   w   h ])



    %%
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

function moveyLabel(a)
    xh = get(gca,'ylabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(1) = p(1)*a ;        % double the distance,
    % negative values put the label below the axis
    set(xh,'position',p);   % set the new position
end

function y = percent95(x)
    I = x>prctile(x, 2.5)& x<prctile(x, 97.5);
    y = x(I);
end
