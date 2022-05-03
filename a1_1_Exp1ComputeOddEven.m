clear
close all

load ('a0_0_merge.mat','dataAllExp1'   )


%%




for iMonkey = 1:2
tic
    data = dataAllExp1{iMonkey};
    data = data(data.trialValid,:);
    sessionDay = day(data.timestamp,  'dayofyear');

    USessionDay = unique(sessionDay);
    nSessionDay  = numel(USessionDay );
    data.iSessionDay = data.ledIntensity*0;
    for iSessionDay = 1:nSessionDay
        I = sessionDay ==  USessionDay(iSessionDay);
        data(I,:).iSessionDay = zeros(nnz(I),1)+ iSessionDay;
    end

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

    image = data.image_number;
    UImage = unique(image);
    nImage  = numel(UImage );


    
    for iImage = 1:nImage
        iImage
        IImage = data.image_number == UImage(iImage);
        nIImage = nnz(IImage);

        for iSessionEvenOdd =  1:2
            I = mod(data.iSessionDay    ,2) == mod(iSessionEvenOdd,2);
            cData = data(IImage & I,:);
            for iStimSite = 1:nLED

                IStimSite = cData.activeLED == ULED(iStimSite);
                cDataSite = cData(IStimSite,:);


                IHit =  cDataSite.trialStim &   cDataSite.monkeyCorrect ;
                IFA  = ~cData.trialStim     & ~cData.monkeyCorrect;
                dPrime(iMonkey).d(iImage,iStimSite,iSessionEvenOdd)  = dPrimeFun(IHit(cDataSite.trialStim ==1), IFA(~cData.trialStim==1));
    
    
            end
        end
    end


end
save('a1_1_Exp1ComputeOddEven.mat', 'dPrime');

