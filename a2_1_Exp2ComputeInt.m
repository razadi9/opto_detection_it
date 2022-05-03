clear
close all

load ('a0_0_merge.mat','dataAllExp2')


%%

nBt = 1e4;


stimReport(2).data =1;



for iMonkey = 1:2

    tic
    data = dataAllExp2{iMonkey};
    data = data(data.trialValid,:);
    image = data.image_number;
    UImage = unique(image);
    nImage  = numel(UImage );



    ULED = unique(data.activeLED);
    ULED(ULED==0) = [];
    nLED = numel(ULED);
    UInt = unique(data.ledIntensity);
    if iMonkey ==1
        UInt(5) = [];
    end
    nInt = numel(UInt);

    NNN(iMonkey).noStim = sum(data.trialStim == 0,1);
    NNN(iMonkey).stim1 = sum(data.trialStim == 1 & data.activeLED == ULED(1),1);
    NNN(iMonkey).stim2 = sum(data.trialStim == 1 & data.activeLED == ULED(2),1);

    for iImage = 1:nImage
        IImage = data.image_number == UImage(iImage);
        for iStimSite = 1:nLED
            IStimSite = data.activeLED == ULED(iStimSite);
            for iInt = 1:nInt
                IInt = data.ledIntensity ==  UInt(iInt);

                ISt = data.trialStim & IImage & IStimSite & IInt;
                INS = ~data.trialStim & IImage;

                IHit =   data.monkeyCorrect(ISt) ;
                IFA  = ~data.monkeyCorrect(INS);
                dPrime(iMonkey).d(iImage,iStimSite,iInt)  =  dPrimeFun(IHit, IFA);


            end

        end
    end


    %
    %% permutation

    dPrime(iMonkey).dPerm = nan(nImage,2,nInt,nBt);

    for iStimSite = 1:nLED
        iStimSite
        IStimSite = data.activeLED == ULED(iStimSite);
        for iInt = 1:nInt
            IInt = data.ledIntensity ==  UInt(iInt);
            ISt = data.trialStim & IStimSite & IInt;
            NSt = find(ISt);
            
            for iImage = 1:nImage


                IImage = data.image_number(NSt) == UImage(iImage);
                n = numel(NSt);
                for iPerm = 1:nBt
                    IHit =   data.monkeyCorrect(NSt(IImage)) ;
                    dPrime(iMonkey).dPerm(iImage,iStimSite,iInt,iPerm)  = nnz(IHit)/numel(IHit);%dPrimeFun(IHit, IFA);
                    NSt = NSt(randperm(n,n));
                end
            end
        end
    end

    if iMonkey ==1
        r1 = 2;
    else
        r1 = 2;
    end


    dPrime(iMonkey).rPerm = nan(2,nBt);
    for iPerm =1:nBt
        for iStimSite = 1:nLED
            r = corr(squeeze(dPrime(iMonkey).dPerm(:,iStimSite,:,iPerm)),'type','Spearman');
            r2 = nan(100,1);
            ii = 0;
            for iR1=r1:size(r,1)
                for iR2=iR1+1:size(r,2)
                    ii = ii + 1;
                    r2(ii) = r(iR1,iR2);
                    %[iR1,iR2]
                end
            end
            dPrime(iMonkey).rPerm(iStimSite,iPerm) = zFisherInv( mean(zFisher(r2(1:ii))));
        end
    end

    toc

    UIntAll{iMonkey} = UInt;


    for iSite = 1:2
        dPerm = dPrime(iMonkey).rPerm(iSite,:);
        dObs = dPrime(iMonkey).rPerm(iSite,1);
        dPrime(iMonkey).p(iSite) = nnz(dObs < dPerm)/numel(dPerm);
    end
end
dPrime(1).rPerm(1:2,1)
dPrime(2).rPerm(1:2,1)
save('a2_1_Exp2ComputeInt.mat', 'dPrime', 'UIntAll', 'NNN');

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

