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
    [h,L,MX,MED,bw] = violin(y,'edgecolor','k','x',1:40,'facecolor',RzTlBx.colorOrder(1));
xlim([.5, 40.5])
    % h:     figure handle
    % L:     Legend handle
    % MX:    Means of groups
    % MED:   Medians of groups
    % bw:    bandwidth of kernel
    legend off


    plot(1:40, dPrime(iMonkey).d(:,iSite),'k.');

    xticks(10:10:100)
    yticks(-2:2:10)
    ylabel('Performance (d'')','FontSize',7,'FontName','Arial')
    xlabel('Image #'  ,'FontSize',7,'FontName','Arial')

    y = ylim;
    ylim(y + diff(y)/20.*[-1 1])
    title('Detection profile')
    axP=[0.05    0.05   0.9   .95 ];
    h = axP(4)*1/3-.075;
    w = axP(3)*4/5-.025;
    axesMain( [axP(1)    1-(axP(2)+h)   w   h ])




    %% Detection Profile - hist

    hAx2 = axes;
    hold on

    dPerm = permute(std(dPrime(iMonkey).dPerm(:,iSite,:),0,1),[3,2,1]);
    dObs = std(dPrime(iMonkey).d(:,iSite,:),0,1);


    histPerm(dPerm, dObs, true, true)

    ylabel('Performance SD (\sigma_{d''})','FontSize',7,'FontName','Arial')
    title('Permutation Test','FontSize',7,'FontName','Arial')


    p = get(hAx1,'Position');
    h = p(4);
    w = axP(3)*1/5-.025;

    axesMain( [axP(1)+axP(3)*4/5    p(2)   w   h ])

    %% Correlation Between and Within
    clear dPrime
    load('a1_1_Exp1ComputeCorr.mat')
    hAx3 = axes;
    hold on


    [h,L,MX,MED,bw, medianLine] = violin( {dPrime(iMonkey).rWithin1, dPrime(iMonkey).rWithin2,dPrime(iMonkey).rBetween}, 'plotmedian', true, 'x',[1,2,4]);
    legend off



    xticks([1,2,4]);
    xticklabels({'Within Site #1', 'Within Site #2', 'Between Sites'})
    yticks(0:.2:1)
    xlim([.5 , 4.5]);
    y = ylim;
    ylim(y + diff(y)/20.*[-1 1])

    %plot([1,5],(yl(1)-.05)*[1,1],'k-')
    ylabel('Correlation coefficient','FontSize',7,'FontName','Arial')
    title('Correlation of stimulation sites')

    h = axP(4)*.8/3-.05;
    w = axP(3)*3/5-.025;
    p = get(hAx1,'Position');
    axesMain( [axP(1)    p(3)-.05-h  w   h ])


    %% Correlation Between and Within - hist

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



    %% Psychometric 
    load('a2_1_Exp2ComputeInt.mat')
    hAx5 = axes;
    hold on

    plot(UInt*1000, squeeze(dPrime(iMonkey).d(:,2,:) )','-','Marker','.','MarkerSize',8)
    xlim([-5,105])
    xlabel('Illumination power (\muW)','FontSize',7,'FontName','Arial')
    ylabel('Performance (d'')','FontSize',7,'FontName','Arial')
    %title('Permutation Test','FontSize',7,'FontName','Arial')
    yticks(-2:2:10)
    xticks(0:50:1000)
    title('Psychometric functions')

    h = axP(4)*1.2/3-.05;
    w = axP(3)*3.5/5-.025;
    p = get(hAx3,'Position');
    axesMain( [axP(1)    axP(1)  w   h ])


    %%
    %% Correlation Between and Within - hist

    hAx6 = axes;
    hold on

    dPerm = dPrime(iMonkey).rPerm(1,:);
    dObs = dPrime(iMonkey).r(1);


    histPerm(dPerm, dObs, true, false)

    ylabel('Average correlation coefficient','FontSize',7,'FontName','Arial')
    title('Permutation Test','FontSize',7,'FontName','Arial')

    y = ylim;
    ylim(y+diff(y)*.1.*[0 1])
    p = get(hAx5,'Position');
    h = p(4);
    w = axP(3)*1.5/5-.025;

    axesMain( [axP(1)+axP(3)*3.5/5    p(2)   w   h ])




    saveas(hFig, 'fig3_Exp1_2___a1_1.pdf')
    saveas(hFig, 'fig3_Exp1_2___a1_1.png')
fakdjfa
end




function histPerm(dPerm, dObs, dObsTextDwon, eText)

    h = histogram(dPerm,'orientation','horizontal','EdgeColor','none','FaceColor',RzTlBx.colorOrder(2),'FaceAlpha',.5);
    set(get(h,'Parent'),'xdir','r')
    set(gca, 'YAxisLocation','right')

    xticks([]);
    yticks(0:.5:10)
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
    xlabel('Frequency','FontSize',7,'FontName','Arial')
    xl = xl+diff(xl)*[0,.05];
    xlim(xl)
    a = get(gca,'XTickLabel');

end
