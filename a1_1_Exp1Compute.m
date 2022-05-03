clear
close all

load ('a0_0_merge.mat','dataAllExp1'   )


%%

nBt = 1e3;
nPrm = 1e4;

stimReport(2).data =1;



for iMonkey = 1:2
    tic
    data = dataAllExp1{iMonkey};
    data = data(data.trialValid,:);
    image = data.image_number;
    UImage = unique(image);
    
    if eGrayExclude
        UImage(UImage==35) = [];
    end

    nImage  = numel(UImage );

    UTrialType = [0,1,2]; %no stim, stim and control



    % merging two control leds for sp since there is not enough data to
    % compute the d primes
    if iMonkey ==2
        nLED = 3;
        I = data.activeLED == 15 ;
        data(I,:).activeLED = ones(nnz(I),1)*18;
    end
    
    ULED = unique(data.activeLED);
    ULED(ULED==0) = [];
    nLED = numel(ULED);


NNN(iMonkey).noStim = sum(data.trialStim == 0 )
if iMonkey ==1
    NNN(iMonkey).catch = sum(data.trialStim == 1 & data.activeLED == ULED(3) );
    NNN(iMonkey).stim1 = sum(data.trialStim == 1 & data.activeLED == ULED(1) );
    NNN(iMonkey).stim2 = sum(data.trialStim == 1 & data.activeLED == ULED(2) );
    I = ~(data.trialStim == 1 & data.activeLED == ULED(3));
    nnz(data(I,:).monkeyCorrect)/nnz(I)

    I1 = data.trialStim == 1 & (data.activeLED == ULED(1) | data.activeLED == ULED(2));
    I2 = ~data.trialStim;
    I = I1 | I2;
    [tableStim,chi2Stim,pStim] = crosstab(data(I,:).trialStim,data(I,:).stimReport)
    sum(tableStim,'all')

    I1 = data.trialStim == 1 & (data.activeLED == ULED(3)                             );
    I2 = ~data.trialStim;
    I = I1 | I2;
    [tableCatch,chi2Catch,pCatch] = crosstab(data(I,:).trialStim,data(I,:).stimReport)
    sum(tableCatch,'all')
    1==1
elseif true
    NNN(iMonkey).catch = sum(data.trialStim == 1 & data.activeLED == ULED(1) );
    NNN(iMonkey).stim1 = sum(data.trialStim == 1 & data.activeLED == ULED(2) );
    NNN(iMonkey).stim2 = sum(data.trialStim == 1 & data.activeLED == ULED(3) );

    I = ~(data.trialStim == 1 & data.activeLED == ULED(1));
    nnz(data(I,:).monkeyCorrect)/nnz(I)


    I1 = data.trialStim == 1 & (data.activeLED == ULED(2) | data.activeLED == ULED(3));
    I2 = ~data.trialStim;
    I = I1 | I2;
    [tableStim,chi2Stim,pStim] = crosstab(data(I,:).trialStim,data(I,:).stimReport)
    sum(tableStim,'all')

    I1 = data.trialStim == 1 & (data.activeLED == ULED(1)                             );
    I2 = ~data.trialStim;
    I = I1 | I2;
    [tableCatch,chi2Catch,pCatch] = crosstab(data(I,:).trialStim,data(I,:).stimReport)
    sum(tableCatch,'all')
    
end

NNN(iMonkey).noStim +NNN(iMonkey).catch + NNN(iMonkey).stim1 +  NNN(iMonkey).stim2
    
    % bootstrap for calculating CIs
    for iImage = 1:nImage
        iImage
        IImage = data.image_number == UImage(iImage);
        cData = data(IImage,:);
        nIImage = nnz(IImage);

        for iStimSite = 1:nLED

            IStimSite = cData.activeLED == ULED(iStimSite);
            cDataSite = cData(IStimSite,:);
            IHit =  cDataSite.trialStim &   cDataSite.monkeyCorrect ;
            IFA  = ~cData.trialStim     & ~cData.monkeyCorrect;
            dPrime(iMonkey).d(iImage,iStimSite)  = dPrimeFun(IHit(cDataSite.trialStim ==1), IFA(~cData.trialStim==1));


            n = size(cData,1 );
            cData0 = cData;
            for iBt = 1:nBt
                N = randi(n,n,1);

                IStimSite = cData(N,:).activeLED == ULED(iStimSite);
                cDataSite = cData(N(IStimSite),:);
                IHit =  cDataSite.trialStim &   cDataSite.monkeyCorrect ;
                IFA  = ~cData(N,:).trialStim & ~cData(N,:).monkeyCorrect;
                
                dPrime(iMonkey).Bt(iImage,iStimSite,iBt)  = dPrimeFun(IHit(cDataSite.trialStim==1), IFA( ~cData(N,:).trialStim ==1));
            end
            dPrime(iMonkey).CI(iImage,iStimSite,1:2) = prctile(dPrime(iMonkey).Bt(iImage,iStimSite,:),[2.5,97.5]);
        end
    end

    %% permutation
    data1 =  data;
    clear img0 IStimSite0
    %
    for iStimSite = 1:nLED+1
        if iStimSite<=nLED
            IStimSite0{iStimSite} = data.activeLED == ULED(iStimSite);
        else
            IStimSite0{iStimSite} = data.activeLED == 0;
        end
        img0{iStimSite} = data( IStimSite0{iStimSite},:).image_number;
    end

    dPrime(iMonkey).dPerm = nan(nImage,nLED, nPrm);

    for iPerm = 1:nPrm
        iPerm
        for iStimSite = 1:nLED+1
            data1( IStimSite0{iStimSite},:).image_number = datasample(img0{iStimSite}, numel(img0{iStimSite}),'Replace',false);
        end


        for iImage = 1:nImage

            IImage = data1.image_number == UImage(iImage);
            cData = data1(IImage,:);
            nIImage = nnz(IImage);

            for iStimSite = 1:nLED

                IStimSite = cData.activeLED == ULED(iStimSite);
                cDataSite = cData(IStimSite,:);
                IHit =  cDataSite.trialStim &   cDataSite.monkeyCorrect ;
                IFA  = ~cData.trialStim & ~cData.monkeyCorrect;
                dPrime(iMonkey).dPerm(iImage,iStimSite, iPerm)  = dPrimeFun(IHit(cDataSite.trialStim==1), IFA(~cData.trialStim ==1));
            end
        end

    end




%{

    %% permutation
    data1 =  data;
    clear stimSite0 IImage0
    %
    for iImage = 1:nImage
        IImage0{iImage} = data1.image_number == UImage(iImage);
            stimSite0{iImage} = data( IImage0{iImage},:).activeLED;

    end



    dPrime(iMonkey).dPermBetweenSites = nan(nImage,nLED, nBt);

    for iPerm = 1:nBt
        iPerm
        for iImage = 1:nImage
            data1( IImage0{iImage},:).activeLED = datasample( stimSite0{iImage}, numel(stimSite0{iImage}),'Replace',false);
        end

        for iStimSite   = 1:2

            IStimSite = data1.activeLED  == ULED(iStimSite);
            cData = data1(IStimSite,:);
            nIStimSite = nnz(IStimSite);

            for iImage = 1:nImage

                IImage = cData.image_number == UImage(iImage);
                cDataImage = cData(IImage,:);
                IHit =  cDataSite.trialStim &   cDataSite.monkeyCorrect ;
                IFA  = ~cData.trialStim & ~cData.monkeyCorrect;
                dPrime(iMonkey).dPermBetweenSites(iImage,iStimSite, iPerm)  = dPrimeFun(IHit(cDataSite.trialStim ==1), IFA(~cData.trialStim==1));
            end
        end        

    end
%}
toc
end
save('a1_1_Exp1Compute.mat', 'dPrime','NNN');

