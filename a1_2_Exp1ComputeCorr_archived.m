clear
close all

load ('a0_0_merge.mat','dataAllExp1'   )


%%

nBt = 1e4;


stimReport(2).data =1;



for iMonkey = 1:2
    tic

    data = dataAllExp1{iMonkey};
    data = data(data.trialValid,:);
    image = data.image_number;
    UImage = unique(image);
    nImage  = numel(UImage );

    UTrialType = [0,1,2]; %no stim, stim and control

    ULED = unique(data.activeLED);
    
    if iMonkey ==1
        ULED([1,4]) = [];
    else
        ULED(1:2) = [];

    end

    nLED = numel(ULED);

    for iStimSite = 1:2
        IStimSite = data.activeLED == ULED(iStimSite);
        for iImage = 1:nImage
            IImage = data.image_number == UImage(iImage);
            NStim{iStimSite,iImage} = find(IStimSite & IImage & data.trialStim);
            NNoStim{iImage}         = find(IImage & ~data.trialStim);
        end
    end

    dPrime(iMonkey).dBt1        = nan(nImage, nBt, 2);
    dPrime(iMonkey).dBt2        = nan(nImage, nBt, 2);
    dPrime(iMonkey).rWithin1    = nan(nBt, 1);
    dPrime(iMonkey).rWithin2    = nan(nBt, 1);
    dPrime(iMonkey).rBetween1    = nan(nBt, 1);
    dPrime(iMonkey).rBetween2    = nan(nBt, 1);
    
    for iImage = 1:nImage
        %iImage

        NStimPrm = [NStim{1,iImage};NStim{2,iImage}];
        nStimPrm(1) = numel(NStim{1,iImage});
        nStimPrm(2) = numel(NStim{2,iImage});

        NNoStimPrm = [NNoStim{iImage};NNoStim{iImage}];
        nNoStimPrm = numel(NNoStim{iImage});
        
        nStimIStimSiteI = nan(nLED, 1);
        nStimIStimSiteI2 = nStimIStimSiteI;
        for iStimSite = 1:nLED
            nStimIStimSiteI(iStimSite) = numel(NStim{iStimSite,iImage});
            nStimIStimSiteI2(iStimSite) = round(nStimIStimSiteI(iStimSite)/2);
        end
        nNNoStim = numel(NNoStim{iImage});
        nNNoStim2 = round(nNNoStim/2);

        % half split
        for iBt = 1:nBt

            for iStimSite = 1:nLED
                n = nStimIStimSiteI(iStimSite);
                n2 = nStimIStimSiteI2(iStimSite);
                N = NStim{iStimSite,iImage}(randperm(n,n));
                IHit1 =  data(N(1:n2),:).trialStim &   data(N(1:n2),:).monkeyCorrect ;
                IHit2 =  data(N(n2+1:end),:).trialStim &   data(N(n2+1:n),:).monkeyCorrect ;

                %{
                n = nNNoStim;
                n2 = nNNoStim2;
                N = NNoStim{iImage} (randperm(n,n));
                IFA1    =  ~data(N(1:n2),:).trialStim &   ~data(N(1:n2),:).monkeyCorrect ;
                IFA2    =  ~data(N(n2+1:end),:).trialStim &   ~data(N(n2+1:end),:).monkeyCorrect ;
                %}

                dPrime(iMonkey).dBt1(iImage, iBt, iStimSite)  = nnz(IHit1)/numel(IHit1);%dPrimeFun(IHit,  IFA);
                dPrime(iMonkey).dBt2(iImage, iBt, iStimSite)  = nnz(IHit2)/numel(IHit2);%dPrimeFun(IHit,  IFA);
            end
        end
    end
    
    for iBt = 1:nBt
        dPrime(iMonkey).rWithin1(iBt) = corr(dPrime(iMonkey).dBt1(:,iBt,1), dPrime(iMonkey).dBt2(:,iBt,1),'type','Spearman');
        dPrime(iMonkey).rWithin2(iBt) = corr(dPrime(iMonkey).dBt1(:,iBt,2), dPrime(iMonkey).dBt2(:,iBt,2),'type','Spearman');
        dPrime(iMonkey).rBetween1(iBt) = corr(dPrime(iMonkey).dBt1(:,iBt,1), dPrime(iMonkey).dBt2(:,iBt,2),'type','Spearman');
        dPrime(iMonkey).rBetween2(iBt) = corr(dPrime(iMonkey).dBt1(:,iBt,2), dPrime(iMonkey).dBt2(:,iBt,1),'type','Spearman');
    end

    
    toc
end

save('a1_2_Exp1ComputeCorr.mat', 'dPrime');



