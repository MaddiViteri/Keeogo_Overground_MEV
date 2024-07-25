function [all_data] = coactivationratio_09_01_rms(all_data,emg_trial_names,muscle_names,forceStruct,trial_name,thresholds)
%% extension vs flexion muscle groups
% thresholds = 50;
emg_freq = 2000;
time = forceStruct.time;
ForceLT = forceStruct.f1(:,2);
ForceRT = forceStruct.f2(:,2);

name = emg_trial_names;
emg_time_full = all_data.(name).time;
%% Data Check
% emg_freq = 2000;
% 
% time = trial.neutral_control05.force.time;
% ForceLT = trial.neutral_control05.force.f1(:,2);
% ForceRT = trial.neutral_control05.force.f2(:,2);
% % 
% % ForceLT = sqrt(sum(ForceLT.^2,2));
% % ForceRT = sqrt(sum(ForceRT.^2,2)); 
% name = emg_trial_names{1};
% emg_time_full = all_data.(name).time;
% thresholds = 50;

%% Analysis
% Match the EMG data with the Trial Data
    if contains(lower(name),lower(trial_name))
        m = 1:length(muscle_names);
        rquad = all_data.(name).rms.(muscle_names{6});
        rham = all_data.(name).rms.(muscle_names{4});
        rgas = all_data.(name).rms.(muscle_names{2});
        lquad = all_data.(name).rms.(muscle_names{5});
        lham = all_data.(name).rms.(muscle_names{3});
        lgas = all_data.(name).rms.(muscle_names{1});

        % % Comparison from EMG set up from other co-contraction- RMS
        % rquad = all_data.(name).rms.(muscle_names{14});
        % rham = all_data.(name).rms.(muscle_names{6});
        % rgas = all_data.(name).rms.(muscle_names{4});
        % lquad = all_data.(name).rms.(muscle_names{13});
        % lham = all_data.(name).rms.(muscle_names{5});
        % lgas = all_data.(name).rms.(muscle_names{3});

        % % Comparison from EMG set up from other co-contraction- Norm EMG
        % rquad = all_data.(name).emg.(muscle_names{14});
        % rham = all_data.(name).emg.(muscle_names{6});
        % rgas = all_data.(name).emg.(muscle_names{4});
        % lquad = all_data.(name).emg.(muscle_names{13});
        % lham = all_data.(name).emg.(muscle_names{5});
        % lgas = all_data.(name).emg.(muscle_names{3});
        
        % % EMG set up from other co-contraction studies- RMS
        % llatgas = all_data.(name).rms.(muscle_names{1});
        % rlatgas = all_data.(name).rms.(muscle_names{2});
        % lmedgas = all_data.(name).rms.(muscle_names{3});
        % rmedgas = all_data.(name).rms.(muscle_names{4});
        % llatham = all_data.(name).rms.(muscle_names{5});
        % rlatham = all_data.(name).rms.(muscle_names{6});
        % lmedham = all_data.(name).rms.(muscle_names{7});
        % rmedham = all_data.(name).rms.(muscle_names{8});
        % llatquad = all_data.(name).rms.(muscle_names{9});
        % rlatquad = all_data.(name).rms.(muscle_names{10});
        % lmedquad = all_data.(name).rms.(muscle_names{11});
        % rmedquad = all_data.(name).rms.(muscle_names{12});
        % lrecfem = all_data.(name).rms.(muscle_names{13});
        % rrecfem = all_data.(name).rms.(muscle_names{14});

        % Combine Flexion/Extension groups
        right_flexion = mean([rham,rgas],2);
        left_flexion = mean([lham,lgas],2);
        right_extension = rquad;
        left_extension = lquad;

        % Combine Flexion/Extension groups
        % right_flexion = mean([rlatham,rmedham,rlatgas,rmedgas],2);
        % left_flexion = mean([llatham,lmedham,llatgas,lmedgas],2);
        % right_extension = mean([rmedquad,rlatquad,rrecfem],2);
        % left_extension = mean([lmedquad,llatquad,lrecfem],2);
        % rgas = mean([rlatgas, rmedgas],2);
        % lgas = mean([llatgas, lmedgas],2);
        % rham = mean([rmedham, rlatham],2);
        % lham = mean([lmedham, llatham],2);

        sum_emg_right = sum([right_flexion,right_extension],2);
        sum_emg_left = sum([left_flexion,left_extension],2);
%% Right Leg
        % Determine which Muscle group is higher for calculation
        for  i = 1:length(right_extension)
            if abs(right_extension(i))>=abs(right_flexion(i))
                all_data.(name).CCI_r(i) = right_flexion(i)*100/(right_extension(i)*100)*(right_flexion(i)*100+right_extension(i)*100);
                all_data.(name).CCR_zero_r(i) = 1 - abs(right_flexion(i))/abs(right_extension(i));
                low(i) = abs(right_flexion(i));
                high(i) = abs(right_extension(i));
            else
                all_data.(name).CCI_r(i) = right_extension(i)*100/(right_flexion(i)*100)*(right_flexion(i)*100+right_extension(i)*100);
                all_data.(name).CCR_zero_r(i) = abs(right_extension(i))/abs(right_flexion(i)) - 1;
                low(i) = abs(right_extension(i));
                high(i) = abs(right_flexion(i));
            end
        end

    % Determine each right step
    [~, heel_strike_data,toe_off_data] = detect_cycles(ForceRT,time,thresholds);
    % Determine minimum number of complete cycles
    nn_r = min(numel(heel_strike_data(:,1)),numel(toe_off_data(:,1)));
    heel_time = cell2mat(heel_strike_data(:,1));
    toe_time = cell2mat(toe_off_data(:,1));
    % Determine the stance time of each step
    for time_index = 1:nn_r
        if heel_time(time_index) > toe_time(time_index)
            for ci = 2:nn_r
                stance_time(ci) = round(1000*(toe_time(ci))) - round(1000*heel_time(ci-1));
            end
        else
            for ci = 1:nn_r
                stance_time(ci) = round(1000*(toe_time(ci))) - round(1000*heel_time(ci));
            end
        end
    end

    % Calculate co-contraction values
    for time_index = 1:nn_r
        % Heel strike index higher than toe off index for same index #
        if heel_time(time_index) > toe_time(time_index)
            for ci = 2:nn_r
                emg_time = heel_time(ci-1):(1/2000):toe_time(ci);
                emg_stance_range = round(round(heel_time(ci-1),3)*emg_freq):round(round(toe_time(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % Calculate co-contraction for each phase
                    all_data.(name).CCI_r_stride_rms(ci) =mean(all_data.(name).CCI_r(emg_stance_range(1:di)));
                    all_data.(name).CCI_r_loading_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCI_r_early_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCI_r_mid_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCI_r_late_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.6):di)));
                    % Summing all the muscles for each stance phase
                    all_data.(name).sum_r_loading_rms(ci) = mean(sum_emg_right(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_r_early_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_r_mid_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_r_late_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.6):di)));
                    % all_data.(name).CCI_r_weight_accept_rms(ci) = trapz(emg_time(1:round(di/5)),all_data.(name).CCI_r(emg_stance_range(1:round(di/5))));
                    % all_data.(name).CCI_r_weight_accept_shi_rms(ci) = trapz(all_data.(name).CCI_r(emg_stance_range(1:round(di/5))))/round(di/5);
                end
                new_t1(1:100,ci) = linspace(heel_time(ci-1),toe_time(ci),100);
                all_data.(name).rms_stance.rgas(:,ci) = interp1(emg_time_full,rgas,new_t1(:,ci)');
                all_data.(name).rms_stance.rham(:,ci) = interp1(emg_time_full,rham,new_t1(:,ci)');
                all_data.(name).rms_stance.rquad(:,ci) = interp1(emg_time_full,right_extension,new_t1(:,ci)');
                % all_data.(name).rms_stance.rgas(:,ci) = interp1(emg_time_full,mean([rlatgas,rmedgas],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.rham(:,ci) = interp1(emg_time_full,mean([rlatham,rmedham],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.rquad(:,ci) = interp1(emg_time_full,mean([rlatquad,rmedquad,rrecfem],2),new_t1(:,ci)');
            end
        else
            % Toe off index higher than heel strike index for same index #
            for ci = 1:nn_r
                emg_time = heel_time(ci):1/2000:toe_time(ci);
                emg_stance_range = round(round(heel_time(ci),3)*emg_freq):round(round(toe_time(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % Calculate co-contraction for each phase
                    all_data.(name).CCI_r_stride_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(1:di)));
                    all_data.(name).CCI_r_loading_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCI_r_early_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCI_r_mid_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCI_r_late_rms(ci) = mean(all_data.(name).CCI_r(emg_stance_range(round(di*0.6):di)));
                    % Summing all the muscles for each stance phase
                    all_data.(name).sum_r_loading_rms(ci) = mean(sum_emg_right(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_r_early_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_r_mid_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_r_late_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.6):di)));
                end
                new_t1(1:100,ci) = linspace(heel_time(ci),toe_time(ci),100);
                all_data.(name).rms_stance.rgas(:,ci) = interp1(emg_time_full,rgas,new_t1(:,ci)');
                all_data.(name).rms_stance.rham(:,ci) = interp1(emg_time_full,rham,new_t1(:,ci)');
                all_data.(name).rms_stance.rquad(:,ci) = interp1(emg_time_full,right_extension,new_t1(:,ci)');
                % all_data.(name).rms_stance.rgas(:,ci) = interp1(emg_time_full,mean([rlatgas,rmedgas],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.rham(:,ci) = interp1(emg_time_full,mean([rlatham,rmedham],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.rquad(:,ci) = interp1(emg_time_full,mean([rlatquad,rmedquad,rrecfem],2),new_t1(:,ci)');
            end
        end
    end
    % mean values for each phase
    all_data.(name).CCI_r_stride_mean_rms = mean(all_data.(name).CCI_r_stride_rms(2:end));
    all_data.(name).CCI_r_loading_mean_rms = mean(all_data.(name).CCI_r_loading_rms(2:end));
    all_data.(name).CCI_r_early_mean_rms = mean(all_data.(name).CCI_r_early_rms(2:end));
    all_data.(name).CCI_r_mid_mean_rms = mean(all_data.(name).CCI_r_mid_rms(2:end));
    all_data.(name).CCI_r_late_mean_rms = mean(all_data.(name).CCI_r_late_rms(2:end));

    all_data.(name).sum_r_mean_loading_rms= mean(all_data.(name).sum_r_loading_rms(2:end));
    all_data.(name).sum_r_mean_early_rms = mean(all_data.(name).sum_r_early_rms(2:end));
    all_data.(name).sum_r_mean_mid_rms = mean(all_data.(name).sum_r_mid_rms(2:end));
    all_data.(name).sum_r_mean_late_rms =mean(all_data.(name).sum_r_late_rms(2:end));
    % std values for each phase
    all_data.(name).CCI_r_stride_std_rms = std(all_data.(name).CCI_r_stride_rms(2:end));
    all_data.(name).CCI_r_loading_std_rms = std(all_data.(name).CCI_r_loading_rms(2:end));
    all_data.(name).CCI_r_early_std_rms = std(all_data.(name).CCI_r_early_rms(2:end));
    all_data.(name).CCI_r_mid_std_rms = std(all_data.(name).CCI_r_mid_rms(2:end));
    all_data.(name).CCI_r_late_std_rms = std(all_data.(name).CCI_r_late_rms(2:end));

    all_data.(name).sum_r_std_loading_rms= std(all_data.(name).sum_r_loading_rms(2:end));
    all_data.(name).sum_r_std_early_rms= std(all_data.(name).sum_r_early_rms(2:end));
    all_data.(name).sum_r_std_mid_rms = std(all_data.(name).sum_r_mid_rms(2:end));
    all_data.(name).sum_r_std_late_rms =std(all_data.(name).sum_r_late_rms(2:end));
    %% Left Leg
    for  i = 1:length(left_extension)
        if abs(left_extension(i))>=abs(left_flexion(i))
            all_data.(name).CCI_l(i) = left_flexion(i)/left_extension(i)*(left_flexion(i)+left_extension(i))*100;
            all_data.(name).CCR_zero_l(i) = 1 - abs(left_flexion(i))/abs(left_extension(i));
            low_l(i) = abs(left_flexion(i));
            high_l(i) = abs(left_extension(i));
        else
            all_data.(name).CCI_l(i) = left_extension(i)/left_flexion(i)*(left_flexion(i)+left_extension(i))*100;
            all_data.(name).CCR_zero_l(i) = abs(left_extension(i))/abs(left_flexion(i)) - 1;
            low_l(i) = abs(left_extension(i));
            high_l(i) = abs(left_flexion(i));
        end
    end
    [~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds);
    nn_l = min(numel(heel_strike_data_l(:,1)),numel(toe_off_data_l(:,1)));
    heel_time_l = cell2mat(heel_strike_data_l(:,1));
    toe_time_l = cell2mat(toe_off_data_l(:,1));
    for time_index = 1:nn_l
        if heel_time_l(time_index) > toe_time_l(time_index)
            for ci = 2:nn_l
                stance_time_l(ci) = round(1000*(toe_time_l(ci))) - round(1000*heel_time_l(ci-1));
            end
        else
            for ci = 1:nn_l
                stance_time_l(ci) = round(1000*(toe_time_l(ci))) - round(1000*heel_time_l(ci));
            end
        end
    end
    %index rms stance field
    % all_data.(name).rms_stance.lquad = NaN(max(stance_time_l)+1,length(nn_l)+1);
    % all_data.(name).rms_stance.lham = NaN(max(stance_time_l)+1,length(nn_l)+1);
    % all_data.(name).rms_stance.lgas = NaN(max(stance_time_l)+1,length(nn_l)+1);
    % all_data.(name).rms_stance.time_l = NaN(max(stance_time_l)+1,length(nn_l)+1);
    for time_index = 1:nn_l
        if heel_time_l(time_index) > toe_time_l(time_index)
            for ci = 2:nn_l
                emg_time = heel_time_l(ci-1):1/2000:toe_time_l(ci);
                emg_stance_range = round(round(heel_time_l(ci-1),3)*emg_freq)+60:round(round(toe_time_l(ci),3)*emg_freq)+60; %includes the 30msec electro mechanical delay
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % Calculate co-contraction for each phase
                    all_data.(name).CCI_l_stride_rms(ci) =mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_loading_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCI_l_early_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCI_l_mid_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCI_l_late_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.6):di)));
                    all_data.(name).CCI_l_mean_stride(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    % Summing all the muscles for each stance phase
                    all_data.(name).sum_l_loading_rms(ci) = mean(sum_emg_left(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_l_early_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_l_mid_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_l_late_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.6):di)));
                end
                new_t1(1:100,ci) = linspace(heel_time_l(ci-1),toe_time_l(ci),100);
                all_data.(name).rms_stance.lgas(:,ci) = interp1(emg_time_full,lgas,new_t1(:,ci)');
                all_data.(name).rms_stance.lham(:,ci) = interp1(emg_time_full,lham,new_t1(:,ci)');
                all_data.(name).rms_stance.lquad(:,ci) = interp1(emg_time_full,left_extension,new_t1(:,ci)');
                % all_data.(name).rms_stance.lgas(:,ci) = interp1(emg_time_full,mean([llatgas,lmedgas],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.lham(:,ci) = interp1(emg_time_full,mean([llatham,lmedham],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.lquad(:,ci) = interp1(emg_time_full,mean([llatquad,lmedquad,lrecfem],2),new_t1(:,ci)');
            end
        else
            for ci = 1:nn_l
                emg_time = heel_time_l(ci):1/2000:toe_time_l(ci);
                emg_stance_range = round(round(heel_time_l(ci),3)*emg_freq):round(round(toe_time_l(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % Calculate co-contraction for each phase
                    all_data.(name).CCI_l_stride_rms(ci) =mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_loading_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCI_l_early_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCI_l_mid_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCI_l_late_rms(ci) = mean(all_data.(name).CCI_l(emg_stance_range(round(di*0.6):di)));
                    all_data.(name).CCI_l_mean_stride(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    % Summing all the muscles for each stance phase
                    all_data.(name).sum_l_loading_rms(ci) = mean(sum_emg_left(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_l_early_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_l_mid_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_l_late_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.6):di)));
                end
                new_t1(1:100,ci) = linspace(heel_time_l(ci),toe_time_l(ci),100);
                all_data.(name).rms_stance.lgas(:,ci) = interp1(emg_time_full,lgas,new_t1(:,ci)');
                all_data.(name).rms_stance.lham(:,ci) = interp1(emg_time_full,lham,new_t1(:,ci)');
                all_data.(name).rms_stance.lquad(:,ci) = interp1(emg_time_full,left_extension,new_t1(:,ci)');
                % all_data.(name).rms_stance.lgas(:,ci) = interp1(emg_time_full,mean([llatgas,lmedgas],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.lham(:,ci) = interp1(emg_time_full,mean([llatham,lmedham],2),new_t1(:,ci)');
                % all_data.(name).rms_stance.lquad(:,ci) = interp1(emg_time_full,mean([llatquad,lmedquad,lrecfem],2),new_t1(:,ci)');
            end
        end
    end

    % all_data.(name).CCI_l_weight_accept_mean_rms = mean(all_data.(name).CCI_l_weight_accept_rms(2:end));
    % all_data.(name).CCR_zero_mean_l_rms = mean(all_data.(name).CCR_zero_l_stride_rms(2:end));
    % Mean For all steps
    all_data.(name).CCI_l_stride_mean_rms = mean(all_data.(name).CCI_l_stride_rms(2:end));
    all_data.(name).CCI_l_loading_mean_rms = mean(all_data.(name).CCI_l_loading_rms(2:end));
    all_data.(name).CCI_l_early_mean_rms = mean(all_data.(name).CCI_l_early_rms(2:end));
    all_data.(name).CCI_l_mid_mean_rms = mean(all_data.(name).CCI_l_mid_rms(2:end));
    all_data.(name).CCI_l_late_mean_rms = mean(all_data.(name).CCI_l_late_rms(2:end));
    % sum
    all_data.(name).sum_l_mean_loading_rms= mean(all_data.(name).sum_l_loading_rms(2:end));
    all_data.(name).sum_l_mean_early_rms = mean(all_data.(name).sum_l_early_rms(2:end));
    all_data.(name).sum_l_mean_mid_rms = mean(all_data.(name).sum_l_mid_rms(2:end));
    all_data.(name).sum_l_mean_late_rms =mean(all_data.(name).sum_l_late_rms(2:end));
    %std
    all_data.(name).CCI_l_stride_std_rms = std(all_data.(name).CCI_l_stride_rms(2:end));
    all_data.(name).CCI_l_loading_std_rms = std(all_data.(name).CCI_l_loading_rms(2:end));
    all_data.(name).CCI_l_early_std_rms = std(all_data.(name).CCI_l_early_rms(2:end));
    all_data.(name).CCI_l_mid_std_rms = std(all_data.(name).CCI_l_mid_rms(2:end));
    all_data.(name).CCI_l_late_std_rms = std(all_data.(name).CCI_l_late_rms(2:end));

    all_data.(name).sum_l_std_loading_rms= std(all_data.(name).sum_l_loading_rms(2:end));
    all_data.(name).sum_l_std_early_rms= std(all_data.(name).sum_l_early_rms(2:end));
    all_data.(name).sum_l_std_mid_rms = std(all_data.(name).sum_l_mid_rms(2:end));
    all_data.(name).sum_l_std_late_rms =std(all_data.(name).sum_l_late_rms(2:end));

    else
        none = 1;
    end
%%
end
