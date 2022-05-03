clear
close all

load ('a0_0_merge.mat','dataAllTraining');


%%

nBoot = 1e4;


stimReport(2).data =[];



for iMonkey = 1:2
    
    data = dataAllTraining{iMonkey};
    sessionDay = day(data.timestamp,  'dayofyear');

    USessionDay = unique(sessionDay);
    nSessionDay  = numel(USessionDay );

    UTrialType = [0,1,2]; %no stim, stim and control

     
    N_trialTypeDay{iMonkey}{ nSessionDay + 1,3} = [];

    for iSessionDay = 1:nSessionDay + 1
        iSessionDay
        for iTrialType = 1:3
            if iSessionDay<nSessionDay + 1
                 N_trialTypeDay{iMonkey}{ iSessionDay, iTrialType} = find(data.trialType == UTrialType(iTrialType) & data.trialValid & sessionDay == USessionDay(iSessionDay));
            else
                % this is for calcuating the last 3 sessions
                 N_trialTypeDay{iMonkey}{ iSessionDay, iTrialType} = find(data.trialType == UTrialType(iTrialType) & data.trialValid & sessionDay >= USessionDay(end-2));
            end



            if ~isempty( N_trialTypeDay{iMonkey}{ iSessionDay, iTrialType}) && numel( N_trialTypeDay{iMonkey}{ iSessionDay, iTrialType}) > 1
                % here we calculate bootstrapped mean, se, and 95% CI of
                % the stim_reports. Then the next mfile will use this data
                % to make the final plots. to keep the results consistent,
                % this bootstarpped data stay saved and will never change.
                stimReport(iMonkey).data{iTrialType, iSessionDay}           = table2array(data( N_trialTypeDay{iMonkey}{ iSessionDay, iTrialType}, 'stimReport'));
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
    y1 = stimReport(iMonkey).data{1, iSessionDay};
    y2 = stimReport(iMonkey).data{2, iSessionDay};
    stimReport(iMonkey).chiSq(iSessionDay,1) = chiSq(y1, y2);
            
    y1 = stimReport(iMonkey).data{1, iSessionDay};
    y2 = stimReport(iMonkey).data{3, iSessionDay};
    stimReport(iMonkey).chiSq(iSessionDay,2) = chiSq(y1, y2);
    end
end
save('a0_1_trainingPhaseCompute2.mat', 'stimReport', 'N_trialTypeDay');



function ch2 = chiSq(y1, y2)

[ch2.tbl,ch2.chi2,ch2.p,ch2.labels] = crosstab([y1;y2],[zeros(size(y1));ones(size(y2))]);


end

