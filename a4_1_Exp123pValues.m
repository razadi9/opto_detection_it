clear
close all



%% correlation between monkeys
load a1_1_Exp1Compute.mat
clear r p n

[r(1),p(1)] = corr(dPrime(1).d(1:38,1),dPrime(2).d(1:38,1),'type', 'Pearson');
[r(2),p(2)] = corr(dPrime(1).d(1:38,1),dPrime(2).d(1:38,2),'type', 'Pearson');
[r(3),p(3)] = corr(dPrime(1).d(1:38,2),dPrime(2).d(1:38,1),'type', 'Pearson');
[r(4),p(4)] = corr(dPrime(1).d(1:38,2),dPrime(2).d(1:38,2),'type', 'Pearson');
min(p)
[min(r),max(r)]
1==1


%% odd/even
load a1_1_Exp1ComputeOddEven.mat
clear r p n
for iMonkey =1:2
    for iSite = 1:size(dPrime(iMonkey).d,2)
        [r(iMonkey,iSite),p(iMonkey,iSite)] = corr(dPrime(iMonkey).d(:,iSite,1),dPrime(iMonkey).d(:,iSite,2),'type','Spearman');
        n(iMonkey,iSite) = numel(dPrime(iMonkey).d(:,iSite,1));
    end
end





%%
load a1_1_Exp1Compute.mat
%% detection profile and gray vs the rest
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


%% Correlation Between and Within


clear dPrime
load('a1_2_Exp1ComputeCorr.mat')
clear p tbl stats
for iMonkey =1:2

    x1 = [dPrime(iMonkey).rWithin1;dPrime(iMonkey).rWithin2];
    x2 = [dPrime(iMonkey).rBetween1;dPrime(iMonkey).rBetween2];

    x = [x1;x2];
    g = [zeros(size(x1));zeros(size(x2))+1];

    [p{iMonkey},tbl{iMonkey},stats{iMonkey}]  = kruskalwallis(x,g,'off');


end


%% Psychometric
load('a2_1_Exp2Compute.mat')

for iMonkey =1:2
    for iStimSite = 1:2
        y = squeeze(dPrime(iMonkey).d(:,iStimSite,:));
        x = meshgrid(1:size(y,2),1:5);
        scatter(x(:),y(:))
        n(iMonkey, iStimSite) = numel(y)-2;
        [r(iMonkey, iStimSite), p(iMonkey, iStimSite)] =corr(x(:),y(:));
    end
end



%%
clear
load a2_1_Exp2ComputeImage.mat
pVal = nan(2,2);
NNoGray = [5,4];
for iMonkey = 1:2
    nImage = size(dPrime(iMonkey).dPerm,1);
    nInt = size(dPrime(iMonkey).dPerm,3);
    nBt = size(dPrime(iMonkey).dPerm,4);
    for iStimSite = 1:2
        x = dPrime(iMonkey).dPerm(:,iStimSite,:,1);
        x = permute(x,[1,3,4,2]);
        x = reshape(x,nImage*nInt,1);
        d1 = std(x,1,1);

        x = dPrime(iMonkey).dPerm(:,iStimSite,:,2:end);
        x = permute(x,[1,3,4,2]);
        x = reshape(x,nImage*nInt,nBt-1);
        dPrm = std(x,1,1);

        pVal(iMonkey,iStimSite) = nnz(d1<dPrm)/numel(dPrm);
    end



end




pValGray = nan(2,2);
IGray = {[false,false,false,false,true],[false,false,false,true,false]};
for iMonkey = 1:2
    for iStimSite = 1:2
        x = dPrime(iMonkey).dPerm(:,iStimSite,:,1);
        x = permute(x,[1,3,4,2]);
        x = [mean(x(~IGray{iMonkey},:),1),x(IGray{iMonkey},:)];
        d1 = std(x,1,2);

        x = dPrime(iMonkey).dPerm(:,iStimSite,:,2:end);
        x = permute(x,[1,3,4,2]);
        x =cat(2,mean(x(~IGray{iMonkey},:,:),1),x(IGray{iMonkey},:,:));
        dPrm = squeeze(std(x,1,2));

        pValGray(iMonkey,iStimSite) = nnz(d1<dPrm)/numel(dPrm);
    end

end


%%
clear
load a2_1_Exp2ComputeInt.mat

for iMonkey = 1:2
    for iStimSite = 1:2
        [r{iMonkey,iStimSite},p{iMonkey,iStimSite}]= corr(squeeze(dPrime(iMonkey).d(:,iStimSite,1:end )),'type','Spearman');
        [min(r{iMonkey,iStimSite}(:)),max(r{iMonkey,iStimSite}(r{iMonkey,iStimSite}<1))]
        
    end
end

1==1
