clear
close all

load ('a0_0_merge.mat','dataAllExp1'   )


%%

eExcludeGray = false;
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

    nAllTrials(iMonkey) = nnz(data.activeLED == ULED(1)|data.activeLED == ULED(2)| ~data.trialStim);
    

    for iStimSite = 1:2
        IStimSite = data.activeLED == ULED(iStimSite);
        for iImage = 1:nImage
            IImage = data.image_number == UImage(iImage);
            NStim{iStimSite,iImage} = find(IStimSite & IImage & data.trialStim);
            NNoStim{iImage}         = find(IImage & ~data.trialStim);
        end
    end

    dPrime(iMonkey).dBt        = nan(nImage, nBt, 2);
    dPrime(iMonkey).dPrm        = nan(nImage, nBt, 2);
    
    for iImage = 1:nImage
        %iImage

        NStimPrm = [NStim{1,iImage};NStim{2,iImage}];
        nStimPrm(1) = numel(NStim{1,iImage});
        nStimPrm(2) = numel(NStim{2,iImage});

        NNoStimPrm = [NNoStim{iImage};NNoStim{iImage}];
        nNoStimPrm = numel(NNoStim{iImage});
        
        for iBt = 1:nBt


            % bootstrap
            for iStimSite = 1:nLED
                n = numel(NStim{iStimSite,iImage});
                N = NStim{iStimSite,iImage} (randi(n,n,1));
                IHit =  data(N,:).trialStim &   data(N,:).monkeyCorrect ;

                n = numel(NNoStim{iImage});
                N = NNoStim{iImage} (randi(n,n,1));
                IFA    =  ~data(N,:).trialStim &   ~data(N,:).monkeyCorrect ;

                dPrime(iMonkey).dBt(iImage, iBt, iStimSite)  = nnz(IHit)/numel(IHit);%dPrimeFun(IHit,  IFA);
            end


            % permuation 


            
            N1 = NStimPrm(1:nStimPrm(1));
            N2 = NStimPrm((1:nStimPrm(2))+nStimPrm(1));
            IHit1 =  data(N1,:).trialStim &   data(N1,:).monkeyCorrect ;
            IHit2 =  data(N2,:).trialStim &   data(N2,:).monkeyCorrect ;

            %{
            n = round(numel(NNoStimPrm)/2);
            N1 = NNoStimPrm(1:n);
            N2 = NNoStimPrm(n+1:end);
            IFA1    =  ~data(N1,:).trialStim &   ~data(N1,:).monkeyCorrect ;
            IFA2    =  ~data(N2,:).trialStim &   ~data(N2,:).monkeyCorrect ;
            %}

            N = NNoStimPrm;
            IFA    =  ~data(N,:).trialStim &   ~data(N,:).monkeyCorrect ;

            dPrime(iMonkey).dPrm(iImage, iBt, 1)  = nnz(IHit1)/numel(IHit1);%dPrimeFun(IHit1,  IFA); %
            dPrime(iMonkey).dPrm(iImage, iBt, 2)  = nnz(IHit2)/numel(IHit2);%dPrimeFun(IHit2,  IFA); %
            
            NStimPrm   = NStimPrm(  randperm(nStimPrm(1)+nStimPrm(2)));
            NNoStimPrm = NNoStimPrm(randperm(nStimPrm(1)+nStimPrm(2)));
        end
    end

    dPrime(iMonkey).rWithin1 = corrlationReport(dPrime(iMonkey).dBt( : ,1:nBt  , 1), dPrime(iMonkey).dBt( :,1:nBt , 1));
    dPrime(iMonkey).rWithin2 = corrlationReport(dPrime(iMonkey).dBt( : ,1:nBt  , 2), dPrime(iMonkey).dBt( :,1:nBt , 2));
    dPrime(iMonkey).rBetween = corrlationReport(dPrime(iMonkey).dBt( : ,1:nBt  , 1), dPrime(iMonkey).dBt( :,1:nBt , 2));

    dPrime(iMonkey).rPrm = nan(nBt, 1);
    dPrime(iMonkey).fitPrm = nan(nBt, 2);
    
    N = 1:40;
    if  eExcludeGray
        N(35) = [];
    end
    for iBt =1:nBt
        dPrime(iMonkey).rPrm(iBt) = corr(dPrime(iMonkey).dPrm(N, iBt, 1), dPrime(iMonkey).dPrm(N, iBt, 2), 'type', 'Pearson');
        dPrime(iMonkey).fitPrm(iBt,1:2) = polyfit(dPrime(iMonkey).dPrm(N, iBt, 1), dPrime(iMonkey).dPrm(N, iBt, 2),1);
    end

    nnz(dPrime(iMonkey).rPrm(1)>dPrime(iMonkey).rPrm(2:end  ))/numel(dPrime(iMonkey).rPrm(2:end  ))
    [dPrime(iMonkey).rPrm(1), median(dPrime(iMonkey).rPrm(2:end  )),nnz(dPrime(iMonkey).rPrm(1)>dPrime(iMonkey).rPrm(2:end  ))/numel(dPrime(iMonkey).rPrm(2:end  ))]    
    toc

end

if ~eExcludeGray
    save('a1_2_Exp1ComputeCorr.mat', 'dPrime','nAllTrials', '-v7.3');
else
    save('a1_2_Exp1ComputeCorrExcludeGray.mat', 'dPrime','nAllTrials', '-v7.3');
end



function r = corrlationReport(d1, d2)
    r0 = corr(d1, d2, 'type','Pearson');
    n = size(r0,1);
    r = nan(n*n/2, 1);
    kR = 0;
    for iR = 1:n
        for jR = iR+1:n
            kR = kR +1;
            r(kR) = r0(iR, jR) ;
        end
    end
    r(kR+1:end) = [];
end


