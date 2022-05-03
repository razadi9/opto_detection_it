clear
%close all

load a1_1_Exp1Compute.mat
e2ndHist = true;


%%
hFig = figure(1);
clf;

set(hFig,'Units','centimeters')
set(hFig, 'Position',[ 15.4869   13.2997   18.3   22])
iH = 0;


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




    %% Detection Profile

    if iMonkey ==1
        NSite = 1;
    else
        NSite = [2,3];
    end
    for iSite = NSite

        iH = iH+1;
        hAx{iH} = axes;
        hold on

        clear y
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
        ylim([min(cellfun(@min, y)),max(cellfun(@max, y))]+[-.15,.1]+[0,.1])



        %ylim([-.1,4.5])
        % h:     figure handle
        % L:     Legend handle
        % MX:    Means of groups
        % MED:   Medians of groups
        % bw:    bandwidth of kernel
        legend off


        plot(1:40, dPrime(iMonkey).d(:,iSite),'k.');

        xticks(10:10:100)

        s = {'Monkey Ph', 'Monkey Sp'};
        ss = {'2','1','2'};
        text(  20  ,max(ylim),[s{iMonkey}, ' - stimulation site #', ss{iH}],'FontSize',9,'FontWeight','normal'    ,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','helvetica','Rotation',0, 'BackgroundColor','none')
        yticks(0:fix(max(ylim)):10)





        y = ylim;
        ylim(y + diff(y)/20.*[-1 1])
        y = ylim;
        if y(1)>0
            ylim(y.*[0,1])
        end
        axP=[0.05    0.05   0.9   .95 ];
        if iMonkey ==1 && iSite ==1
            iH=  1;
        elseif iMonkey ==2 && iSite ==2
            iH = 2;
        elseif iMonkey ==2 && iSite ==3
            iH = 3;
        else
            error
        end

        if iH ==1
            h = axP(4)*1/3-.075;
            w = axP(3)*4/5-.025;
            axesMain( [axP(1)    1-(axP(2)+h)   w   h ])
            title('Detection profile')

        else
            p = get(hAx{iH-1}, 'position');
            axesMain( [p(1)    p(2)-p(4)-.075   p(3)   p(4) ])

        end

        if iH ==3
            xlabel('Image index'  ,'FontSize',9,'FontName','Arial')
        end
        if iH ==2
            ylabel('Performance (d'')','FontSize',9,'FontName','Arial')

        end
        x= xlim;
        y= xlim;

        if iMonkey == 1 && iSite == 1
            n = NNN(iMonkey).(['stim',num2str(iSite)])+NNN(iMonkey).noStim;
        elseif iMonkey == 2 && iSite == 2
            n = NNN(iMonkey).(['stim',num2str(1)])+NNN(iMonkey).noStim;
        elseif iMonkey == 2 && iSite == 3
            n = NNN(iMonkey).(['stim',num2str(2)])+NNN(iMonkey).noStim;
        end
        s = num2strComma(n);
        text(min(xlim)+diff(xlim)*.03*0,min(ylim)+diff(ylim)*.05*0, sprintf('  n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

        set(get(gca, 'title'),'fontweight','normal')

        p = get(gca,'Position');

        %axesText = axes('Parent', gcf, ...
        %    'Units', 'normalized', ...
        %    'Position', [0, 0, 1, 1], ...
        %    'Visible', 'off', ...
        %    'XLim', [0, 1], ...
        %    'YLim', [0, 1], ...
        %    'NextPlot', 'add');

        %axes(axesText);
        %text(p(1)-.02,p(2)+p(4)+.02, 'a', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');


        %% Detection Profile - hist

        hAx2 = axes;
        hold on

        dPerm = permute(std(dPrime(iMonkey).dPerm(:,iSite,:),0,1),[3,2,1]);
        dObs = std(dPrime(iMonkey).d(:,iSite,:),0,1);


        histPerm(dPerm, dObs, true, iH==1)
        yticks(.10:.7:10)
        ylim([.1,.8])


        xh = get(gca,'ylabel'); % handle to the label object
        p = get(xh,'position'); % get the current position property
        p(1) = p(1)*1.5 ;        % double the distance,
        % negative values put the label below the axis
        set(xh,'position',p);   % set the new position


        if iH ==2
            ylabel('Performance SD (\sigma_{d''})','FontSize',9,'FontName','Arial')
            ylabel('Performance SD','FontSize',9,'FontName','Arial')
        end
        if iH == 1
            %title('Permutation Test','FontSize',9,'FontName','Arial')
        end

        p = get(hAx{iH},'Position');
        h = p(4);
        w = axP(3)*1/5-.025;

        axesMain( [axP(1)+axP(3)*4/5    p(2)   w   h ])



        p = get(gca,'Position');
        %axes(axesText);
        %text(p(1)-.02,p(2)+p(4)+.02, 'b', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');

    end
end



saveas(hFig, 'figED1_Exp1___a1_91.fig')
saveas(hFig, 'figED1_Exp1___a1_91.pdf')
saveas(hFig, 'figED1_Exp1___a1_91.png')




%% Correlation Between and Within
clear dPrime
load('a1_2_Exp1ComputeCorr.mat')
dPrime(iMonkey).rWithin1 = percent95(dPrime(iMonkey).rWithin1);
dPrime(iMonkey).rWithin2 = percent95(dPrime(iMonkey).rWithin2);
dPrime(iMonkey).rBetween = percent95(dPrime(iMonkey).rBetween);

%%

clear hAx
hFig = figure(2);
clf;
set(hFig,'Units','centimeters')
set(hFig, 'Position',[ 15.4869   13.2997   8.8   6.6])

hAx1 = axes;
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
hAx1.XTickLabel = tickLabels;

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

p = get(hAx1,'Position');
axesMain( [p(1)    p(2)  p(3)*1.45/(1.45+.8)-.025   p(4) ])


s = num2strComma(nAllTrials(iMonkey));

text(min(xlim)+diff(xlim)*.03*0,min(ylim)+diff(ylim)*.05*0, sprintf('  n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

moveyLabel(-1.5);
set(gca,'TickLength',get(gca,'TickLength')*3)
p = get(gca,'Position');



p = get(hAx1, 'Position');
a  = axes;
rgbImage = imread('array_schematic_beveled_SP.png','BackgroundColor',[1 1 1]);
rgbImage = fliplr(rgbImage);
rgbImage = imresize(rgbImage,.5);
%rgbImage = fliplr(rgbImage);
imshow(rgbImage);




set(gca,'position',[p(1)+p(3)*.05, p(2)+p(4)*(-.05)+.1, .1*22.0133/6.5969.*.8, .1*18.3000/8.7842.*.8])
box off
axis off


1==1

%% Correlation Between and Within - hist
hAx2 = axes;
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


p = get(hAx1,'Position');
h = p(4);
w = axP(3)*(1- (1.45/(1.45+.8))-.025)-.025;

axesMain( [p(1)+p(3)+.025, p(2),   w   h ])


moveyLabel(0.5);









saveas(hFig, 'figED2_Exp1___a1_91.fig')
saveas(hFig, 'figED2_Exp1___a1_91.pdf')
saveas(hFig, 'figED2_Exp1___a1_91.png')

%% Psychometric
close all

hFig = figure(3);

%%

clf;
set(hFig,'Units','centimeters')
set(hFig, 'Position',[ 15.4869   13.2997   18.3/2   22])

iH = 0;
e2ndHist = true;
load('a2_1_Exp2ComputeIntFit.mat')
iH =0;
for iMonkey = 1:2
    if iMonkey == 1
        N = 1;
    else
        N = 1:2;
    end
    for iSite = N



        iH = iH+1;
        hAx5{iH} = axes;
        hold on
        if iH ==1
            h = axP(4)*1/3-.075;
            w = axP(3)*1.65/(1.65+.8)-.025;

            axesMain( [axP(1)+.05    1-(axP(2)+h)   w-0.05   h ])
            


        else
            p = get(hAx5{iH-1}, 'position');
            axesMain( [p(1)    p(2)-p(4)-.075   p(3)   p(4) ])

        end

            title('Intensity psychometric functions')
            set(get(gca, 'title'),'fontweight','normal')
        s = {'Monkey Ph', 'Monkey Sp'};
        ss = {'2','1','2'};
        title(  [s{iMonkey}, ' - stimulation site #', ss{iH}])


        %if e2ndHist
        %    w = axP(3)*2/5-.025*2;
        %else
        %    w = axP(3)*3/5-.025*2;
        %end
        hold on
        if iH == 1
            xlim([-5,105])
        else
            xlim([-5,155])
        end

        %axesMain( [axP(1)+axP(3)*2/5+.025    p(2)   w   h ])
        p = get(gca,'Position');
        axesMain( p)

        if iMonkey ==1
            N1 = 1:4;
            N2 = 5;
        else
            N1 = [1,2,3,5];
            N2 = 4;
        end
        h = plot(UIntAll{iMonkey}*1000, squeeze(dPrime(iMonkey).d(iSite,N1,:) )','-','Marker','o','MarkerSize',5,'LineWidth',1);
        for iHH= 1:numel(h)
            set(h(iHH),'color',RzTlBx.colorOrder(iHH+1))
        end
        plot(UIntAll{iMonkey}*1000, squeeze(dPrime(iMonkey).d(iSite, N2,:) )','-','Marker','s','MarkerSize',5,'LineWidth',1,'color',colNoImage)
        x = UIntAll{iMonkey}(4:5)*1000;
        y= squeeze(dPrime(iMonkey).d(iSite, N2,4:5) );

        a  = atand((diff(y)./diff(x))./(diff(ylim)./diff(xlim)))-10;
        text(mean(x), mean(y)-.2,{'no image'},'Color',colNoImage, 'BackgroundColor', 'none','Rotation', a,'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 7,'FontName', 'helvetica');
        if iH ==3
            xlabel('Illumination power (\muW)','FontSize',9,'FontName','Arial')
        elseif iH == 2
            ylabel('Performance (d'')','FontSize',9,'FontName','Arial')
            moveyLabel(.75);

        end
        %title('Permutation Test','FontSize',7,'FontName','Arial')
        y = fix(ylim);
        yticks([0,y(2)])
        
        %xticks(0:50:1000)

        n = NNN(iMonkey).noStim + NNN(iMonkey).(['stim',num2str(iSite)]);
        y = ylim;
        ylim([-.5,4])
        s = num2strComma(n);
        text(max(xlim)-diff(xlim)*.05*0,min(ylim)+diff(ylim)*.00, sprintf('n = %s trials',s),'FontSize',8,'FontWeight','normal'    ,'HorizontalAlignment','right','VerticalAlignment','bottom','FontName','helvetica','Rotation',0, 'BackgroundColor','none')

        xticks(0:50:1000)
        1==1
        %h = axP(4)*1.2/3-.05;
        %w = axP(3)*3.5/5-.025;
        %p = get(hAx3,'Position');
        %axesMain( [axP(1)    axP(1)  w   h ])

        p = get(gca,'Position');
        %axes(axesText);
        %text(p(1)-.02,p(2)+p(4)+.02, 'd', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize',12, 'FontName','arial','FontWeight','bold');

        %%
        %%  hist
        %
        if e2ndHist
            hAx6 = axes;
            hold on

            histPerm(dPrime(iMonkey).s(iSite,2:end), dPrime(iMonkey).s(iSite,1), false, false)
            y = ylim.*[0,1];
            ylim(y)
            ylim([0,.3])
            y = ylim;
            yticks(y)
            ylabel('Coefficient SD','FontSize',7,'FontName','Arial')
            set(gca,'TickLength',get(gca,'TickLength')*3)

            moveyLabel(.5);


            p = get(hAx5{iH},'Position');
            h = p(4);
            w = axP(3)*.8/(.8+1.65)-.025*3;

            axesMain( [p(1)+p(3)+.025*2, p(2),   w   h ])








        end


    end
    %}
end

saveas(hFig, 'figED3_Exp2___a1_91.fig')
saveas(hFig, 'figED3_Exp2___a1_91.pdf')
saveas(hFig, 'figED3_Exp2___a1_91.png')
fakdjfa







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
