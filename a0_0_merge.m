% this code will read the data files and convert them to new unified
% matrixes for each experiment. these matrixes will include all the data
% for both monkeys.

% note that in the archive_0, monkey 1 is sprocket and monkey 2 is phelyx,
% here I changed the order to be consistent with the paper.

clear
close all

%% ====================== training phase ======================

for iMonkey = 1:2

    % read the data
    if iMonkey == 1       % Phelyx
        load('data/phelyx_training_data_uniform.mat');
        trainingData = renamevars(trainingData,{'broken_fixation_trial', 'correction_loop_trial_nostim_trial' , 'correction_loop_trial_stim_trial', 'dim_dot_opacity', 'dim_dots', 'stim_trial', 'monkey_correct'}, ...
            {'brokenFixation'       , 'correctionLoop_noStim'              , 'correctionLoop_stim'             , 'cueOpacity'     , 'cueTrial', 'trialStim' , 'monkeyCorrect' });

        trainingData.stimReport = ~trainingData.looked_up;

        %INoStim = trainingData.trialType == 0;
        %IStim = trainingData.trialType == 1;
        %IControl = trainingData.trialType == 2;
        % LED number:
        %    37: Catch
        %    13,17: Virus
    else
        
        trainingData =  readtable('data/experiment_1_sprocket_IT_detection_logfile.csv');
        trainingData.Properties.VariableNames = {'timestamp','blocknumber','broken_fixation_trial','stim_trial','side','image_number','looked_up','monkey_correct','juice_amount','rewarded','ledIntensity','stimLength','time_state_system_start','time_trial_start','time_accept_trial','time_of_stimulation','time_mid_stim_temp_read','time_response_targets_show','time_response_made','time_end_trial','temp_at_start_right','temp_at_start_left','temp_before_stim_right','temp_before_stim_left','temp_mid_stim_right','temp_mid_stim_left','temp_after_stim_right','temp_after_stim_left','temp_at_end_of_trial_right','temp_at_end_of_trial_left','temp_after_choice_right','temp_after_choice_left','channel_1_intensity','channel_2_intensity','channel_3_intensity','channel_4_intensity','channel_5_intensity','channel_6_intensity','channel_7_intensity','channel_8_intensity','channel_9_intensity','channel_10_intensity','channel_11_intensity','channel_12_intensity','channel_13_intensity','channel_14_intensity','channel_15_intensity','channel_16_intensity','channel_17_intensity','channel_18_intensity','channel_19_intensity','channel_20_intensity','channel_21_intensity','channel_22_intensity','channel_23_intensity','channel_24_intensity','channel_33_intensity','channel_34_intensity','channel_35_intensity','channel_36_intensity','channel_37_intensity','channel_38_intensity','channel_39_intensity','channel_40_intensity','channel_41_intensity','channel_42_intensity','channel_43_intensity','channel_44_intensity','channel_45_intensity','channel_46_intensity','channel_47_intensity','channel_48_intensity','channel_49_intensity','channel_50_intensity','channel_51_intensity','channel_52_intensity','channel_53_intensity','channel_54_intensity','channel_55_intensity','channel_56_intensity','timed_out_trial','eyetracker_sentinal','looked_up_correction_loop_accum','looked_down_correction_loop_accum','downward_correction_loop_trial','upward_correction_loop_trial'};

        trainingData = renamevars(trainingData,{'broken_fixation_trial','looked_up_correction_loop_accum','looked_down_correction_loop_accum','downward_correction_loop_trial','upward_correction_loop_trial', 'monkey_correct' ,  'stim_trial'}, ...
                                       {'brokenFixation'       , 'correctionLoop1'               , 'correctionLoop2'                 , 'correctionLoop3'              ,  'correctionLoop4'           , 'monkeyCorrect'  ,  'trialStim' });


        trainingData.correctionLoop_stim  = false(size(trainingData,1),1);
        trainingData.correctionLoop_noStim  = false(size(trainingData,1),1);

        I = ~cellfun(@isempty,trainingData.correctionLoop1) | cellfun(@isempty,trainingData.correctionLoop2) | cellfun(@isempty,trainingData.correctionLoop3) | cellfun(@isempty,trainingData.correctionLoop4);
        trainingData(I,:).correctionLoop_stim = cellfun(@(x) ~strcmp(x,'0'), trainingData(I,:).correctionLoop1) & cellfun(@(x) ~strcmp(x,'0'), trainingData(I,:).correctionLoop2);
        trainingData(I,:).correctionLoop_noStim = cellfun(@(x) ~strcmp(x,'0'), trainingData(I,:).correctionLoop3) & cellfun(@(x) ~strcmp(x,'0'), trainingData(I,:).correctionLoop4);

        trainingData.stimReport = ~trainingData.looked_up;


        trainingData.trialType = trainingData.trialStim;
        I = trainingData.trialType  == 1;

        trainingData(I,:).trialType = trainingData(I,:).side+1;
    end

    trainingData.trialValid = trainingData.brokenFixation == 0 & trainingData.correctionLoop_stim == 0 & trainingData.correctionLoop_noStim == 0;
    dataAllTraining{iMonkey} = trainingData;
end

save('a0_0_merge.mat','dataAllTraining')

clear trainingData dataAllTraining
%% ====================== Experiment 1 ======================
dataAllExp1{2} = [];


for iMonkey =1:2
    if iMonkey == 1
        load('data/Exp1/allData.mat');
        exp1Data = allData;
        clear allData
        exp1Data = renamevars(exp1Data,{'broken_fixation_trial', 'correction_loop_trial_nostim_trial' , 'correction_loop_trial_stim_trial', 'dim_dot_opacity', 'dim_dots', 'stim_trial', 'monkey_correct'}, ...
                                           {'brokenFixation'   , 'correctionLoop_noStim'              , 'correctionLoop_stim'             , 'cueOpacity'     , 'cueTrial', 'trialStim' , 'monkeyCorrect' });
                
    else
        exp1Data =  readtable('data/Exp1/sprocket_40img_datafile.csv');
        exp1Data.Properties.VariableNames = {'timestamp', 'blocknumber', 'broken_fixation_trial', 'stim_trial', 'side','image_number', 'looked_up', 'monkey_correct','juice_amount', 'rewarded', 'ledIntensity', 'stimLength', 'time_state_system_start', 'time_trial_start', 'time_accept_trial', 'time_of_stimulation', 'time_mid_stim_temp_read', 'time_response_targets_show', 'time_response_made', 'time_end_trial', 'temp_at_start_right', 'temp_at_start_left', 'temp_before_stim_right', 'temp_before_stim_left', 'temp_mid_stir_right', 'temp_mid_stim_left', 'temp_after_stim_right', 'temp_after_stim_left', 'temp_at_end_of_trial_right', 'temp_at_end_of_trial_left', 'temp_after_choice_right', 'temp_after_choice_left', 'channel_1_intensity', 'channel_2_intensity', 'channel_3_intensity', 'channel_4_intensity', 'channel_5_intensity', 'channel_6_intensity', 'channel_7_intensity', 'channel_8_intensity','channel_9_intensity', 'channel_10_intensity', 'channel_11_intensity', 'channel_12_intensity', 'channel_13_intensity', 'channel_14_intensity', 'channel_15_intensity', 'channel_16_intensity', 'channel_17_intensity', 'channel_18_intensity', 'channel_19_intensity', 'channel_20_intensity', 'channel_21_intensity', 'channel_22_intensity', 'channel_23_intensity', 'channel_24_intensity', 'channel_33_intensity', 'channel_34_intensity', 'channel_35_intensity', 'channel_36_intensity', 'channel_37_intensity', 'channel_38_intensity', 'channel_39_intensity', 'channel_40_intensity', 'channel_41_intensity', 'channel_42_intensity', 'channel_43_intensity', 'channel_44_intensity', 'channel_45_intensity', 'channel_46_intensity', 'channel_47_intensity', 'channel_48_intensity', 'channel_49_intensity', 'channel_50_intensity', 'channel_51_intensity', 'channel_52_intensity', 'channel_53_intensity', 'channel_54_intensity', 'channel_55_intensity', 'channel_56_intensity', 'timed_out_trial', 'eyetracker_sentinel', 'looked_up_correction_accum', 'looked_down_correction_accum', 'downward_correction_trial', 'upward_correction_trial', 'state_number', 'active_LED'};

        exp1Data = renamevars(exp1Data,{'broken_fixation_trial', 'upward_correction_trial' , 'downward_correction_trial',  'stim_trial', 'monkey_correct', 'active_LED'}, ...
                                           {'brokenFixation'   , 'correctionLoop_noStim'    , 'correctionLoop_stim'     , 'trialStim' , 'monkeyCorrect'  , 'activeLED' });
    end


    exp1Data.stimReport = exp1Data.looked_up;
    exp1Data.trialValid = exp1Data.brokenFixation == 0 & exp1Data.correctionLoop_stim == 0 & exp1Data.correctionLoop_noStim == 0;
    dataAllExp1{iMonkey} = exp1Data;

end
save('a0_0_merge.mat','dataAllExp1','-append')
clear dataAllExp1 exp1Data






%% ====================== Experiment 2 ======================


if false
if iMonkey == 1
    allData =  readtable('experiment_1_sprocket_IT_detection_logfile_psychometric_functions.csv');
    allData.Properties.VariableNames = {'timestamp', 'blocknumber', 'broken_fixation_trial', 'stim_trial', 'side','image_number', 'looked_up', 'monkey_correct','juice_amount', 'rewarded', 'ledIntensity', 'stimLength', 'time_state_system_start', 'time_trial_start', 'time_accept_trial', 'time_of_stimulation', 'time_mid_stim_temp_read', 'time_response_targets_show', 'time_response_made', 'time_end_trial', 'temp_at_start_right', 'temp_at_start_left', 'temp_before_stim_right', 'temp_before_stim_left', 'temp_mid_stir_right', 'temp_mid_stim_left', 'temp_after_stim_right', 'temp_after_stim_left', 'temp_at_end_of_trial_right', 'temp_at_end_of_trial_left', 'temp_after_choice_right', 'temp_after_choice_left', 'channel_1_intensity', 'channel_2_intensity', 'channel_3_intensity', 'channel_4_intensity', 'channel_5_intensity', 'channel_6_intensity', 'channel_7_intensity', 'channel_8_intensity','channel_9_intensity', 'channel_10_intensity', 'channel_11_intensity', 'channel_12_intensity', 'channel_13_intensity', 'channel_14_intensity', 'channel_15_intensity', 'channel_16_intensity', 'channel_17_intensity', 'channel_18_intensity', 'channel_19_intensity', 'channel_20_intensity', 'channel_21_intensity', 'channel_22_intensity', 'channel_23_intensity', 'channel_24_intensity', 'channel_33_intensity', 'channel_34_intensity', 'channel_35_intensity', 'channel_36_intensity', 'channel_37_intensity', 'channel_38_intensity', 'channel_39_intensity', 'channel_40_intensity', 'channel_41_intensity', 'channel_42_intensity', 'channel_43_intensity', 'channel_44_intensity', 'channel_45_intensity', 'channel_46_intensity', 'channel_47_intensity', 'channel_48_intensity', 'channel_49_intensity', 'channel_50_intensity', 'channel_51_intensity', 'channel_52_intensity', 'channel_53_intensity', 'channel_54_intensity', 'channel_55_intensity', 'channel_56_intensity', 'timed_out_trial', 'eyetracker_sentinel', 'looked_up_correction_accum', 'looked_down_correction_accum', 'downward_correction_trial', 'upward_correction_trial', 'state_number', 'active_LED'};
else
    load allData2.mat
    I = allData.image_number == 4;
    
end
end



dataAllExp2{2} = [];


for iMonkey =1:2
    if iMonkey == 1
        load('data/Exp2/allData2.mat');
        exp2Data = allData;
        clear allData
        exp2Data = renamevars(exp2Data,{'broken_fixation_trial', 'correction_loop_trial_nostim_trial' , 'correction_loop_trial_stim_trial', 'dim_dot_opacity', 'dim_dots', 'stim_trial', 'monkey_correct'}, ...
                                           {'brokenFixation'       , 'correctionLoop_noStim'              , 'correctionLoop_stim'             , 'cueOpacity'     , 'cueTrial', 'trialStim' , 'monkeyCorrect' });
                
    else
    exp2Data =  readtable('data/Exp2/experiment_1_sprocket_IT_detection_logfile_psychometric_functions.csv');
    exp2Data.Properties.VariableNames = {'timestamp', 'blocknumber', 'broken_fixation_trial', 'stim_trial', 'side','image_number', 'looked_up', 'monkey_correct','juice_amount', 'rewarded', 'ledIntensity', 'stimLength', 'time_state_system_start', 'time_trial_start', 'time_accept_trial', 'time_of_stimulation', 'time_mid_stim_temp_read', 'time_response_targets_show', 'time_response_made', 'time_end_trial', 'temp_at_start_right', 'temp_at_start_left', 'temp_before_stim_right', 'temp_before_stim_left', 'temp_mid_stir_right', 'temp_mid_stim_left', 'temp_after_stim_right', 'temp_after_stim_left', 'temp_at_end_of_trial_right', 'temp_at_end_of_trial_left', 'temp_after_choice_right', 'temp_after_choice_left', 'channel_1_intensity', 'channel_2_intensity', 'channel_3_intensity', 'channel_4_intensity', 'channel_5_intensity', 'channel_6_intensity', 'channel_7_intensity', 'channel_8_intensity','channel_9_intensity', 'channel_10_intensity', 'channel_11_intensity', 'channel_12_intensity', 'channel_13_intensity', 'channel_14_intensity', 'channel_15_intensity', 'channel_16_intensity', 'channel_17_intensity', 'channel_18_intensity', 'channel_19_intensity', 'channel_20_intensity', 'channel_21_intensity', 'channel_22_intensity', 'channel_23_intensity', 'channel_24_intensity', 'channel_33_intensity', 'channel_34_intensity', 'channel_35_intensity', 'channel_36_intensity', 'channel_37_intensity', 'channel_38_intensity', 'channel_39_intensity', 'channel_40_intensity', 'channel_41_intensity', 'channel_42_intensity', 'channel_43_intensity', 'channel_44_intensity', 'channel_45_intensity', 'channel_46_intensity', 'channel_47_intensity', 'channel_48_intensity', 'channel_49_intensity', 'channel_50_intensity', 'channel_51_intensity', 'channel_52_intensity', 'channel_53_intensity', 'channel_54_intensity', 'channel_55_intensity', 'channel_56_intensity', 'timed_out_trial', 'eyetracker_sentinel', 'looked_up_correction_accum', 'looked_down_correction_accum', 'downward_correction_trial', 'upward_correction_trial', 'state_number', 'active_LED'};

        exp2Data = renamevars(exp2Data,{'broken_fixation_trial', 'upward_correction_trial' , 'downward_correction_trial',  'stim_trial', 'monkey_correct', 'active_LED'}, ...
                                           {'brokenFixation'    , 'correctionLoop_noStim'         , 'correctionLoop_stim'     , 'trialStim' , 'monkeyCorrect' ,'activeLED'});
    end
    exp2Data.trialValid = exp2Data.brokenFixation == 0 & exp2Data.correctionLoop_stim == 0 & exp2Data.correctionLoop_noStim == 0;
    dataAllExp2{iMonkey} = exp2Data;

end
save('a0_0_merge.mat','dataAllExp2','-append')
clear dataAllExp2 exp2Data
%% ====================== Experiment 3 ======================
dataAllExp3{2} = [];
for iMonkey =1:2
    if iMonkey == 1
        load('data/Exp3/DataE3M2');
    else
        load('data/Exp3/DataE3M1');
    end
    allData.trialValid  = allData.broken_fixation_trial == 0 & allData.correction_loop_trial_nostim_trial==0 & allData.correction_loop_trial_stim_trial ==0 ;

    dataAllExp3{iMonkey} = allData;
end

save('a0_0_merge.mat','dataAllExp3','-append')

