clear
close all

load ('a0_0_merge.mat','dataAllExp2')


%%

myfittype = fittype('a*x.^b','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
options = fitoptions(myfittype);
options.Lower = [1e-9,1e-9];
options.Upper = [Inf,1];
options.StartPoint = [1,.5];

nPerm = 1e4;
eGrayExclude = false;


stimReport(2).data =1;


dPrime(2).s = nan(2,nPerm);
tic
for iMonkey = 1:2

    tic
    data = dataAllExp2{iMonkey};
    data = data(data.trialValid,:);
    image = data.image_number;
    UImage = unique(image);


    if eGrayExclude
        if iMonkey ==1
            I = data.image_number == UImage(5);
        else
            I = data.image_number == UImage(4);
        end
        data(I,:) = [];
        image = data.image_number;
        UImage = unique(image);
        

    end
    
    nImage  = numel(UImage );
    nData = size(data,1);


    ULED = unique(data.activeLED);
    ULED(ULED==0) = [];
    nLED = numel(ULED);
    UInt = unique(data.ledIntensity);
    if iMonkey ==1
        UInt(5) = [];
    end
    nInt = numel(UInt);

    x = sqrt(UInt);

    for iPerm = 1:nPerm
        for iStimSite = 1:nLED
            IStimSite = data.activeLED == ULED(iStimSite);
            for iImage = 1:nImage
                IImage = data.image_number == UImage(iImage);
                for iInt = 1:nInt
                    IInt = data.ledIntensity ==  UInt(iInt);

                    ISt = data.trialStim & IImage & IStimSite & IInt;
                    INS = ~data.trialStim & IImage;

                    IHit =  data.monkeyCorrect(ISt);
                    IFA  = ~data.monkeyCorrect(INS);
                    d(iImage, iInt)   =  dPrimeFun(IHit, IFA);
                end
                %[fitobject{iImage}, gof(iImage)] = fit(UInt, d(iImage,:)', myfittype, options);
                %y = d(iImage,:)';
                %fitobject{iImage}.a = x\y;
                %gof(iImage).rsquare = 1-sum((x*fitobject{iImage}.a-y).^2)/sum((y-mean(y)).^2);
            end
            y = d(1:nImage,:)';
            [fitobject, gof] = fit(repmat(UInt,nImage,1), y(:), myfittype, options);
            dPrime(iMonkey).fita( iStimSite, iPerm    ) = fitobject.a;
            dPrime(iMonkey).fitb( iStimSite, iPerm    ) = fitobject.b;
            x = (fitobject.a*UInt.^fitobject.b);

            for iImage = 1:nImage

            
            

                y = d(iImage,:)';
                dPrime(iMonkey).a( iStimSite, iPerm,iImage) = x\y;
                dPrime(iMonkey).r2(iStimSite, iPerm,iImage) = 1-sum((x* dPrime(iMonkey).a(iStimSite, iPerm,iImage)-y).^2)/sum((y-mean(y)).^2);

            end
            if iPerm == 1
                dPrime(iMonkey).d(iStimSite, 1:nImage,1:nInt) = d;
            end

            N = data.image_number(IStimSite);
            data.image_number(IStimSite) = N(randperm(numel(N)));
        end
        N = data.image_number(~data.trialStim);
        data.image_number(~data.trialStim ) = N(randperm(numel(N)));
    end
    UIntAll{iMonkey} = UInt;

    toc
end

for iMonkey = 1:2
    for iStimSite = 1:2
        s = std(dPrime(iMonkey).a(iStimSite,:,:),0,3);
        dPrime(iMonkey).s(iStimSite,:) = s;
        (nnz(s(1)<s))/nPerm
    end
end
if ~eGrayExclude
    save('a2_1_Exp2ComputeIntFit.mat', 'dPrime','UIntAll');
else
    save('a2_1_Exp2ComputeIntFitGrayExclude.mat', 'dPrime','UIntAll');
end


