clear
close all
load ('a0_0_merge.mat','dataAllExp3');

nBt = 1e4;
% 1 : time stamp
% 2 : blocknumber
% 3 : broken_fixation_trial
% 4 : correction_loop_trial_stim_trial
% 5 : correction_loop_trial_nostim_trial
% 6 : stim_trial
% 7 : dim_dots
% 8 : dim_dot_opacity
% 9 : array_num
% 10 : image_number
% 11 : looked_up
% 12 : monkey_correct
% 13 : juice_amount
% 15 : ledIntensity
% 27 : temp_before_stim_array1
% 31 : temp_after_stim_array1
% 36 - 84 : channel intensities
% 87 : trialType
% 88 : activeLED



for iMonkey =2:2
    data = dataAllExp3{iMonkey};
    data = data(data.trialValid,:);
    numel(unique(day(data.timestamp,  'dayofyear')))
    size(data,1)
    nnz(data.monkey_correct)/size(data,1)
    
    NNN(iMonkey) = size(data,1);
    UImage = unique(data.image_number);
    nImage = numel(UImage);
    dPrime(iMonkey).Bt = nan(nImage, nBt);

    for iImage = 1:nImage


        IImage = data.image_number == UImage(iImage) ;
        cData = data(IImage,:);

        IHit =  cData.stim_trial &  cData.monkey_correct;
        IFA  = ~cData.stim_trial & ~cData.monkey_correct;

        dPrime(iMonkey).d(iImage)  = dPrimeFun(IHit(cData.stim_trial==1), IFA(~cData.stim_trial==1));

        n = size(cData,1 );
        for iBt = 1:nBt
            N = randi(n,n,1);
            IHit =  cData(N,:).stim_trial &  cData(N,:).monkey_correct;
            IFA  = ~cData(N,:).stim_trial & ~cData(N,:).monkey_correct;
            dPrime(iMonkey).Bt(iImage,iBt) = dPrimeFun(IHit(cData(N,:).stim_trial==1), IFA(~cData(N,:).stim_trial==1));
        end
        dPrime(iMonkey).CIBt(iImage,1:2) = prctile(dPrime(iMonkey).Bt(iImage,1:nBt),[2.5,97.5]);
    end

if iMonkey ==1
    img(iMonkey).index = [4:-1:1;9:-1:6;13:-1:10;17:-1:14;21:-1:18];
    img(iMonkey).g = 5;
elseif iMonkey ==2


    x = [1 5 10 14 18]';x =fliplr([x,x+1,x+2,x+3]);
    img(iMonkey).index = x;
    img(iMonkey).g = 9;
end    

end

save('a3_1_Exp3Compute.mat', 'dPrime','img', 'NNN');

%{



IStim = allData.stim_trial;
ICorr = allData.monkey_correct;

IHit  = IStim  & ICorr;
ICR   = ~IStim & ICorr;
IFA   = ~IStim & ~ICorr;
IMiss = IStim & ~ICorr;




%% dprime profile (observed)
nBt = 1e5;

dPrime     = nan(nImage, nBt);


ISt =  IValid &  IStim;
INS =  IValid & ~IStim;

NSt =  find(ISt);
NNS =  find(INS);

for iImage = 1:nImage
    iImage
    IImage = allData.image_number == UImage(iImage);
    N1 = find(ISt&IImage);
    N2 = find(INS&IImage);
    n1 = nnz(N1);
    n2 = nnz(N2);
    for iBt = 1:nBt
        if iBt==1

            dPrime(iImage,iBt) = d_prime(IHit(N1), IFA(N2), true);
        else

            dPrime(iImage,iBt) = d_prime(IHit(datasample( N1, n1, 'Replace', true)), IFA(datasample( N2, n2, 'Replace', true)), false);
        end
    end
end


%%

if iMonkey ==1
    load a4_1_Exp3_Gray_M1.mat
    x = [1 5 10 14 18]';x =fliplr([x,x+1,x+2,x+3]);
    xG = 9;
elseif iMonkey ==2

    %load a4_1_Exp3_Gray_M2
    load a4_1_Exp3_Gray_M2.mat
    x = [4:-1:1;9:-1:6;13:-1:10;17:-1:14;21:-1:18];
    xG = 5;

end


clf
cla reset

hold on
xlabel('Visual Simulus Intesity')
ylabel('Performance (d'')')
title('Monkey #2')
RzTlBx.boldFig;

dPrimeAll = 0;
for iGrayness = 1:size(x,1)
    N = x(iGrayness,:);
    dPrimeAll =  dPrimeAll+dPrime(N,:);

    %dPrimeImage(

end


dPrimeAll = dPrimeAll/iGrayness;
dPrimeAll(5,:) = dPrime(xG,:);


if nBt >1e3
    n=1e3;
else
    n = nBt;
end
plot(0:4,dPrimeAll([5,1:4],1:n)','Color',RzTlBx.colorOrder(1,.01),'LineWidth',1)


xx = [];
yy = [];
for iGrayness = 1:size(x,1)
    N = x(iGrayness,:);
    plot(1:4,dPrime(N,1),'.-','LineWidth',.5,'Color',RzTlBx.colorOrder(iGrayness+1),'MarkerSize',10)
    xx = [xx;(1:4)'];
    yy = [yy;dPrime(N,1)];
end

[rc,pc] = corr(xx,yy,'type', 'Spearman');

p  =sum(dPrimeAll(4,:)>=dPrimeAll([3,2,1,5],:),2)./size(dPrimeAll,2);
p(p>.5) = 1-p(p>.5);
p = 2*p;
q = mafdr(p,'BH', true);


y = max(yy);

for iP = 1:4
    y = y + .15;
    x = [4-iP,4];
    pValueDraw(x, y*[1,1], q(iP));
end

if pc >=.001
    error
else
    s = sprintf('r = %.2f\np < 0.001***', rc);
end
text(.2, 3, s, 'BackgroundColor', 'none', 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'left', 'FontSize',12);




xticks(0:5)
xlim([-.1,4.1])
%ylim([.8,3])
yticks(0:10);

x = 1:4;
y = dPrimeAll(1:4,1);
yEU = prctile(dPrimeAll(1:4,:)',97.5)'-y;
yEL = y-prctile(dPrimeAll(1:4,:)',2.5)';
errorbar(x,y,yEU, yEL,'o-','Color',RzTlBx.colorOrder(1),'LineWidth',2,'MarkerFaceColor',RzTlBx.colorOrder(1),'MarkerSize',8);
x = 0;
y = dPrimeAll(5,1);
yEU = prctile(dPrimeAll(5,:)',97.5)'-y;
yEL = y-prctile(dPrimeAll(5,:)',2.5)';
errorbar(x,y,yEU, yEL,'o-','Color',RzTlBx.colorOrder(1),'LineWidth',2,'MarkerFaceColor',RzTlBx.colorOrder(1),'MarkerSize',8);


y = ylim;
ylim([y(1), 3.7])


if false
    a4_1_Exp3_Gray_M2
    saveas(gcf,'fig4.fig')
end


jidfalakhjdflajf

%%

function dPrime = d_prime(IHit, IFA, eWarning)
if nargin<3
    eWarning = true;
end

ht = nnz(IHit)/numel(IHit);
fa = nnz(IFA)/numel(IFA);

if ht == 1
    ht = .99;
    if eWarning
        warning('Hit rate = 100%')
    end
elseif ht == 0
    ht = .01;
    if eWarning
        warning('Hit rate = 0%')
    end
end
if fa == 1
    fa = .99;
    if eWarning
        warning('FA rate = 100%')
    end
elseif fa == 0
    fa = .01;
    if eWarning
        warning('FA rate = 0%')
    end
end
dPrime = norminv(ht) - norminv(fa);
end

function [N1,N2] = sampDiv(N)

N = datasample(N, numel(N), 'Replace', false);
n = round(numel(N)/2);
N1 = N(1:n);
N2 = N(n+1:end);
end

function rZ = zFisher(r)
I = r >.99;
r(I) = .99;
I = r <-.99;
r(I) = -.99;
rZ = 0.5*log((1+r)./(1-r));
end

function r = zFisherInv(rZ)
r = (exp(rZ/.5)  - 1)./(1 +  exp(rZ/.5)) ;
end
%}