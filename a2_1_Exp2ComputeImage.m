clear
close all

load ('a0_0_merge.mat','dataAllExp2')


%%

nBt = 1e3;


stimReport(2).data =1;



for iMonkey = 1:2
tic
    data = dataAllExp2{iMonkey};
    data = data(data.trialValid,:);
    numel(unique(day(data.timestamp,  'dayofyear')))
    size(data,1)
    nnz(data.monkeyCorrect)/size(data,1)
    image = data.image_number;
    UImage = unique(image);
    nImage  = numel(UImage );

    UTrialType = [0,1,2]; %no stim, stim and control

    ULED = unique(data.activeLED);
    ULED(ULED==0) = [];
    nLED = numel(ULED);
    UInt = unique(data.ledIntensity);
    if iMonkey ==1
        UInt(5) = [];
    end
    nInt = numel(UInt);


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

    for iStimSite = 1:nLED
        for iInt1 = 1:nInt
            for iInt2 = 1:nInt
                r(iInt1,iInt2) = corr(dPrime(iMonkey).d(:,iStimSite,iInt1),dPrime(iMonkey).d(:,iStimSite,iInt2),'type','Spearman');
            end
        end
        dPrime(iMonkey).r(iStimSite) = zFisherInv( mean(zFisher(r(2:end,2:end)),'all'));
    end
    %
    %% permutation


    for iStimSite = 1:nLED
        iStimSite
        IStimSite = data.activeLED == ULED(iStimSite);
        for iImage = 1:nImage
                IImage = data.image_number ==  UImage(iImage);
            for iInt = 1:nInt




                for iPerm = 1:nBt

                    ISt = data.trialStim & IStimSite& IImage;
                    INS = ~data.trialStim & IImage;

                    if iPerm ==1
                        ISt = find(ISt& data.ledIntensity== UInt(iInt));
                        INS = find(INS);
                    else
                        n = nnz(ISt & data.ledIntensity== UInt(iInt));
                        ISt = datasample(find(ISt),n,'replace',false);
                        n = nnz(INS );
                        INS = datasample(find(INS),n,'replace',false);
                    end


                    IHit =   data.monkeyCorrect(ISt) ;
                    IFA  = ~data.monkeyCorrect(INS);
                    dPrime(iMonkey).dPerm(iImage,iStimSite,iInt,iPerm)  = nnz(IHit)/numel(IHit);%dPrimeFun(IHit, IFA);

                end
            end
        end
    end



    %%



toc


end

save('a2_1_Exp2ComputeImage.mat', 'dPrime');

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

