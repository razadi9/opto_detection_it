clear
close all

load ('a0_0_merge.mat','dataAllTraining');


%%

nBoot = 1e4;


stimReport(2).data =1;



for iMonkey = 1:1
    
    data = dataAllTraining{iMonkey};
    sessionDay = day(data.timestamp,  'dayofyear');

    USessionDay = unique(sessionDay);
    nSessionDay  = numel(USessionDay );

    UTrialType = [0,1,2]; %no stim, stim and control

    clear N_trialTypeDay 

     
    N_trialTypeDay{3} = [];

    for iSessionDay = 1:nSessionDay + 1
        iSessionDay
        for iTrialType = 1:3
            if iSessionDay<nSessionDay + 1
                N_trialTypeDay{iTrialType} = find(data.trialType == UTrialType(iTrialType) & data.trialValid & sessionDay == USessionDay(iSessionDay));
            else
                % this is for calcuating the last 3 sessions
                N_trialTypeDay{iTrialType} = find(data.trialType == UTrialType(iTrialType) & data.trialValid & sessionDay >= USessionDay(end-2));


            end
            if ~isempty(N_trialTypeDay{iTrialType})
                % here we calculate bootstrapped mean, se, and 95% CI of
                % the stim_reports. Then the next mfile will use this data
                % to make the final plots. to keep the results consistent,
                % this bootstarpped data stay saved and will never change.
                stimReport(iMonkey).data{iTrialType, iSessionDay}           = table2array(data(N_trialTypeDay{iTrialType}, 'stimReport'));
                stimReport(iMonkey).btData(iTrialType, iSessionDay,1:nBoot) = bootstrp(nBoot, @mean, stimReport(iMonkey).data{iTrialType, iSessionDay});
                stimReport(iMonkey).m(iTrialType, iSessionDay)              = mean( stimReport(iMonkey).data{iTrialType, iSessionDay});
                stimReport(iMonkey).btSE(iTrialType, iSessionDay)           = std(stimReport(iMonkey).btData(iTrialType, iSessionDay,1:nBoot),1,3);
                stimReport(iMonkey).btCI(iTrialType, iSessionDay,1:2)       = prctile(stimReport(iMonkey).btData(iTrialType, iSessionDay,1:nBoot),[2.5,97.5],3);
            else

                stimReport(iMonkey).data{iTrialType, iSessionDay}           = nan;
                stimReport(iMonkey).btData(iTrialType, iSessionDay,1:nBoot) = nan;
                stimReport(iMonkey).m(iTrialType, iSessionDay)              = nan;
                stimReport(iMonkey).btSE(iTrialType, iSessionDay)           = nan;
                stimReport(iMonkey).btCI(iTrialType, iSessionDay,1:2)       = nan;
                

            end
        end

         
                
        
    


        %permutation test
        for iComp = 2:-1:1
            if iComp == 1
                y1 = stimReport(iMonkey).data{1, iSessionDay};
                y2 = stimReport(iMonkey).data{2, iSessionDay};
            else
                y1 = stimReport(iMonkey).data{1, iSessionDay};
                y2 = stimReport(iMonkey).data{3, iSessionDay};
            end
            y = [y1;y2];
            n = cumsum([numel(y1);numel(y2)]);
            
            d0 = mean(y(1:n(1)))- mean(y(n(1)+1:n(2)));
            d = nan(nBoot, 1);
            for iPrm = 1:nBoot
                y = y(randperm(n(2)));
                d(iPrm) = mean(y(1:n(1)))- mean(y(n(1)+1:n(2)));
            end
            p = (nnz(d0>d) + nnz(d0==d)/2)/nBoot; % the second term is in case all (or most) d are same as d0
            if p>.5; p=1-p;end
            p = p*2;
    
            permTest(iMonkey).d(iSessionDay, iComp,1:nBoot) = d;
            permTest(iMonkey).p(iSessionDay, iComp) = p;
        end
    end
end
save('a0_1_trainingPhaseCompute.mat', 'stimReport','permTest');



%{
            if any(I)
                dPrime(iSessionDay,iTrialType,iBoot)= nnz(stimReport_trialType(I))/numel(stimReport_trialType(I));
                
                 %[dPrime(iSessionDay,iTrialType,iBoot), ht, fa, criterion] = d_prime(IHitB(I), IFAB(I), false);
                 %1==1
                
            else
                dPrime(iSessionDay,iTrialType,iBoot) = nan;
            end
        end
    end
    
end


        n = numel(N);


    for iBoot = 1:nBoot
        if iBoot ==1
            NB = N;
        else
            NB = datasample(N, n,  'replace', true);
        end
        
        sessionDay_trialType =  sessionDay(NB);
        end
        %}
