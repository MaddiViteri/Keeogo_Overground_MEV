function [all_data] = coactivationratio_zero(all_data,emg_trial_names,muscle_names,forceStruct,trial_name,thresholds)
%%

emg_freq = 2000;
time = forceStruct.time;
ForceLT = forceStruct.f1(:,2);
ForceRT = forceStruct.f2(:,2);
name = emg_trial_names;
% thresholds = 50;
% time = trial.neutral_control.force.time;
% ForceLT = trial.neutral_control.force.f1(:,2);
% ForceRT = trial.neutral_control.force.f2(:,2);
% name = emg_trial_names{2};
% trial_name =trial_name{1};
    if contains(lower(name),lower(trial_name))
        m = 1:length(muscle_names);
        rquad = all_data.(name).emg.(muscle_names{6});
        rham = all_data.(name).emg.(muscle_names{4});
        rgas = all_data.(name).emg.(muscle_names{2});
        lquad = all_data.(name).emg.(muscle_names{5});
        lham = all_data.(name).emg.(muscle_names{3});
        lgas = all_data.(name).emg.(muscle_names{1});
        
        % EMG For comparison set up
        % rquad = all_data.(name).rms.(muscle_names{6});
        % rham = all_data.(name).rms.(muscle_names{4});
        % rgas = all_data.(name).rms.(muscle_names{2});
        % lquad = all_data.(name).rms.(muscle_names{5});
        % lham = all_data.(name).rms.(muscle_names{3});
        % lgas = all_data.(name).rms.(muscle_names{1});
        % rquad = all_data.(name).emg.(muscle_names{14});
        % rham = all_data.(name).emg.(muscle_names{6});
        % rgas = all_data.(name).emg.(muscle_names{4});
        % lquad = all_data.(name).emg.(muscle_names{13});
        % lham = all_data.(name).emg.(muscle_names{5});
        % lgas = all_data.(name).emg.(muscle_names{3});
        
        % EMG For comparison set up 2
        % llatgas = all_data.(name).emg.(muscle_names{1});
        % rlatgas = all_data.(name).emg.(muscle_names{2});
        % lmedgas = all_data.(name).emg.(muscle_names{3});
        % rmedgas = all_data.(name).emg.(muscle_names{4});
        % llatham = all_data.(name).emg.(muscle_names{5});
        % rlatham = all_data.(name).emg.(muscle_names{6});
        % lmedham = all_data.(name).emg.(muscle_names{7});
        % rmedham = all_data.(name).emg.(muscle_names{8});
        % llatquad = all_data.(name).emg.(muscle_names{9});
        % rlatquad = all_data.(name).emg.(muscle_names{10});
        % lmedquad = all_data.(name).emg.(muscle_names{11});
        % rmedquad = all_data.(name).emg.(muscle_names{12});
        % lrecfem = all_data.(name).emg.(muscle_names{13});
        % rrecfem = all_data.(name).emg.(muscle_names{14});
        
        right_flexion = mean([rham,rgas],2);
        left_flexion = mean([lham,lgas],2);
        right_extension = rquad;
        left_extension = lquad;

        % % Grouped values for Comparison
        % right_flexion = mean([rlatham,rmedham,rlatgas,rmedgas],2);
        % left_flexion = mean([llatham,lmedham,llatgas,lmedgas],2);
        % right_extension = mean([rmedquad,rlatquad,rrecfem],2);
        % left_extension = mean([lmedquad,llatquad,lrecfem],2);

        sum_emg_right = sum([right_flexion,right_extension],2);
        sum_emg_left = sum([left_flexion,left_extension],2);
        
        %% Right Leg  
        % CCI/R calculation     
        for i = 1:length(right_extension) 
            if abs(right_extension(i))>=abs(right_flexion(i))
                all_data.(name).CCI_r(i) = right_flexion(i)/right_extension(i)*(right_flexion(i)+right_extension(i));
                all_data.(name).CCR_zero_r(i) = 1 - abs(right_flexion(i))/abs(right_extension(i));
                low(i) = abs(right_flexion(i));
                high(i) = abs(right_extension(i));
            else
                all_data.(name).CCI_r(i) = right_extension(i)/right_flexion(i)*(right_flexion(i)+right_extension(i));
                all_data.(name).CCR_zero_r(i) = abs(right_extension(i))/abs(right_flexion(i)) - 1;
                low(i) = abs(right_extension(i));
                high(i) = abs(right_flexion(i));
            end
        end

    %% Right leg stance phase calculation
    [~, heel_strike_data,toe_off_data] = detect_cycles(ForceRT,time,thresholds); % positive (pos) and negative (neg) slope crossing points
    nn_r = min(numel(heel_strike_data(:,1)),numel(toe_off_data(:,1)));
    heel_time = cell2mat(heel_strike_data(:,1));
    toe_time = cell2mat(toe_off_data(:,1));
    for time_index = 1:nn_r
        if heel_time(time_index) > toe_time(time_index)
            for ci = 2:nn_r
                emg_time = heel_time(ci-1):(1/2000):toe_time(ci);
                emg_stance_range = round(round(heel_time(ci-1),3)*emg_freq):round(round(toe_time(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % all_data.(name).CCI_winters_r_stride(ci) =...
                    %     2*100*trapz(emg_time(1:di),low(emg_stance_range(1:di)))/...
                    %     (trapz(emg_time(1:di),low(emg_stance_range(1:di))) + ...
                    %     trapz(emg_time(1:di),high(emg_stance_range(1:di))));
                    % all_data.(name).CCI_r_mean_stride(ci) = mean(all_data.(name).CCI_r(emg_stance_range(1:di)));
                    % all_data.(name).CCI_r_mean_stride(ci) = trapz(emg_time(1:di),all_data.(name).CCI_r(emg_stance_range(1:di)));
                    % all_data.(name).CCI_r_weight_accept(ci) = trapz(emg_time(1:round(di/5)),all_data.(name).CCI_r(emg_stance_range(1:round(di/5))));
                    % all_data.(name).CCI_r_weight_accept_shi(ci) = trapz(all_data.(name).CCI_r(emg_stance_range(1:round(di/5))))/round(di/5);
                    % all_data.(name).CCR_zero_mean_r_stride(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(1:di)));
                    % all_data.(name).CCR_zero_mean_r_stride(ci) =trapz(emg_time(1:di),CCR_zero_r(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_r_stride(ci) =mean(all_data.(name).CCR_zero_r(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_r_loading(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCR_zero_r_early(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCR_zero_r_mid(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCR_zero_r_late(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.6):di)));

                    all_data.(name).sum_r_loading_rms(ci) = mean(sum_emg_right(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_r_early_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_r_mid_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_r_late_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.6):di)));
                end
            end
        else
            for ci = 2:nn_r
                emg_time = heel_time(ci):1/2000:toe_time(ci);
                emg_stance_range = round(round(heel_time(ci),3)*emg_freq):round(round(toe_time(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    % all_data.(name).CCI_winters_r_stride(ci) =...
                    %     2*100*trapz(emg_time(1:di),low(emg_stance_range(1:di)))/...
                    %     (trapz(emg_time(1:di),low(emg_stance_range(1:di))) + ...
                    %     trapz(emg_time(1:di),high(emg_stance_range(1:di))));
                    % all_data.(name).CCI_r_mean_stride(ci) = mean(all_data.(name).CCI_r(emg_stance_range(1:di)));
                    % all_data.(name).CCI_r_mean_stride(ci) = trapz(emg_time(1:di),all_data.(name).CCI_r(emg_stance_range(1:di)));
                    % all_data.(name).CCI_r_weight_accept(ci) = trapz(emg_time(1:round(di/5)),all_data.(name).CCI_r(emg_stance_range(1:round(di/5))));
                    % all_data.(name).CCR_zero_mean_r_stride(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(1:di)));
                    % all_data.(name).CCR_zero_mean_r_stride(ci) =trapz(emg_time(1:di),all_data.(name).CCR_zero_r(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_r_stride(ci) =mean(all_data.(name).CCR_zero_r(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_r_loading(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCR_zero_r_early(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCR_zero_r_mid(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCR_zero_r_late(ci) = mean(all_data.(name).CCR_zero_r(emg_stance_range(round(di*0.6):di)));

                    all_data.(name).sum_r_loading_rms(ci) = mean(sum_emg_right(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_r_early_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_r_mid_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_r_late_rms(ci) = mean(sum_emg_right(emg_stance_range(round(di*0.6):di)));

                end
            end
        end
    end
    %%
    % all_data.(name).CCI_winters_r_mean = mean(all_data.(name).CCI_winters_r_stride(2:end));
    % all_data.(name).CCI_r_mean = mean(all_data.(name).CCI_r_mean_stride(2:end));
    % all_data.(name).CCI_r_weight_accept_mean = mean(all_data.(name).CCI_r_weight_accept(2:end));
    all_data.(name).CCR_zero_mean_r = mean(all_data.(name).CCR_zero_r_stride(2:end));

    all_data.(name).CCR_zero_r_loading_mean = mean(all_data.(name).CCR_zero_r_loading(2:end));
    all_data.(name).CCR_zero_r_early_mean = mean(all_data.(name).CCR_zero_r_early(2:end));
    all_data.(name).CCR_zero_r_mid_mean = mean(all_data.(name).CCR_zero_r_mid(2:end));
    all_data.(name).CCR_zero_r_late_mean = mean(all_data.(name).CCR_zero_r_late(2:end));
    %std
    all_data.(name).CCR_zero_r_loading_std = std(all_data.(name).CCR_zero_r_loading(2:end));
    all_data.(name).CCR_zero_r_early_std = std(all_data.(name).CCR_zero_r_early(2:end));
    all_data.(name).CCR_zero_r_mid_std = std(all_data.(name).CCR_zero_r_mid(2:end));
    all_data.(name).CCR_zero_r_late_std = std(all_data.(name).CCR_zero_r_late(2:end));
    %             new_t1_r(1:100,ci) = linspace(t0_pos1(ci-1),t0_neg1(ci),100)';
    % all_data.(name).CCI = sum(all_data.(name).CCI)*100/length(all_data.(name).CCI);
    % all_data.(name).CCI_winters = 2*100*trapz(time,low)/...
    %     (trapz(time,low) + trapz(time,high));
    %% Left Leg
    % CCI/CCR calculated for entire trial
    for i = 1:length(left_extension)
        if abs(left_extension(i))>=abs(left_flexion(i))
            all_data.(name).CCI_l(i) = left_flexion(i)/left_extension(i)*(left_flexion(i)+left_extension(i));
            all_data.(name).CCR_zero_l(i) = 1 - abs(left_flexion(i))/abs(left_extension(i));
            low_l(i) = abs(left_flexion(i));
            high_l(i) = abs(left_extension(i));
        else
            all_data.(name).CCI_l(i) = left_extension(i)/left_flexion(i)*(left_flexion(i)+left_extension(i));
            all_data.(name).CCR_zero_l(i) = abs(left_extension(i))/abs(left_flexion(i)) - 1;
            low_l(i) = abs(left_extension(i));
            high_l(i) = abs(left_flexion(i));
        end
    end
    % CCI/CCR calculated for stance phase
    [~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds);
    nn_l = min(numel(heel_strike_data_l(:,1)),numel(toe_off_data_l(:,1)));
    heel_time_l = cell2mat(heel_strike_data_l(:,1));
    toe_time_l = cell2mat(toe_off_data_l(:,1));
    for time_index = 1:nn_l
        % if statement depending on which is reported first (heel or toe)
        if heel_time_l(time_index) > toe_time_l(time_index)
            for ci = 2:nn_l
                emg_time = heel_time_l(ci-1):1/2000:toe_time_l(ci);
                emg_stance_range = round(round(heel_time_l(ci-1),3)*emg_freq)+60:round(round(toe_time_l(ci),3)*emg_freq)+60; %includes the 30msec electro mechanical delay
                for di = min(numel(emg_time),numel(emg_stance_range))
                    all_data.(name).CCI_winters_l_stride(ci) =...
                        2*100*trapz(emg_time(1:di),low_l(emg_stance_range(1:di)))/...
                        (trapz(emg_time(1:di),low_l(emg_stance_range(1:di))) + ...
                        trapz(emg_time(1:di),high_l(emg_stance_range(1:di))));
                    % all_data.(name).CCI_l_mean_stride(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_mean_stride(ci) = trapz(emg_time(1:di),all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_weight_accept(ci) = trapz(emg_time(1:round(di/5)),all_data.(name).CCI_l(emg_stance_range(1:round(di/5))));
                    % all_data.(name).CCR_zero_mean_l_stride(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    % all_data.(name).CCR_zero_mean_l_stride(ci) =trapz(emg_time(1:di),all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_l_stride(ci) =mean(all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_l_loading(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCR_zero_l_early(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCR_zero_l_mid(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCR_zero_l_late(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.6):di)));

                    all_data.(name).sum_l_loading_rms(ci) = mean(sum_emg_left(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_l_early_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_l_mid_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_l_late_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.6):di)));
                end
            end
        else
            for ci = 2:nn_l
                emg_time = heel_time_l(ci):1/2000:toe_time_l(ci);
                emg_stance_range = round(round(heel_time_l(ci),3)*emg_freq):round(round(toe_time_l(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_stance_range))
                    all_data.(name).CCI_winters_l_stride(ci) =...
                        2*100*trapz(emg_time(1:di),low_l(emg_stance_range(1:di)))/...
                        (trapz(emg_time(1:di),low_l(emg_stance_range(1:di))) + ...
                        trapz(emg_time(1:di),high_l(emg_stance_range(1:di))));
                    % all_data.(name).CCI_l_mean_stride(ci) = mean(all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_mean_stride(ci) = trapz(emg_time(1:di),all_data.(name).CCI_l(emg_stance_range(1:di)));
                    all_data.(name).CCI_l_weight_accept(ci) = trapz(emg_time(1:round(di/5)),all_data.(name).CCI_l(emg_stance_range(1:round(di/5))));
                    % all_data.(name).CCR_zero_mean_l_stride(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    % all_data.(name).CCR_zero_mean_l_stride(ci) =trapz(emg_time(1:di),all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_l_stride(ci) =mean(all_data.(name).CCR_zero_l(emg_stance_range(1:di)));
                    all_data.(name).CCR_zero_l_loading(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).CCR_zero_l_early(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).CCR_zero_l_mid(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).CCR_zero_l_late(ci) = mean(all_data.(name).CCR_zero_l(emg_stance_range(round(di*0.6):di)));

                    all_data.(name).sum_l_loading_rms(ci) = mean(sum_emg_left(emg_stance_range(1:round(di*0.15))));
                    all_data.(name).sum_l_early_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.15):round(di*0.4))));
                    all_data.(name).sum_l_mid_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.4):round(di*0.6))));
                    all_data.(name).sum_l_late_rms(ci) = mean(sum_emg_left(emg_stance_range(round(di*0.6):di)));
                end
            end
        end
    end
    all_data.(name).CCI_winters_l_mean = mean(all_data.(name).CCI_winters_l_stride(2:end));
    all_data.(name).CCI_l_mean = mean(all_data.(name).CCI_l_mean_stride(2:end));
    all_data.(name).CCI_l_weight_accept_mean = mean(all_data.(name).CCI_l_weight_accept(2:end));
    all_data.(name).CCR_zero_mean_l = mean(all_data.(name).CCR_zero_l_stride(2:end));
    
    %mean
    all_data.(name).CCR_zero_l_loading_mean = mean(all_data.(name).CCR_zero_l_loading(2:end));
    all_data.(name).CCR_zero_l_early_mean = mean(all_data.(name).CCR_zero_l_early(2:end));
    all_data.(name).CCR_zero_l_mid_mean = mean(all_data.(name).CCR_zero_l_mid(2:end));
    all_data.(name).CCR_zero_l_late_mean = mean(all_data.(name).CCR_zero_l_late(2:end));
    %std
    all_data.(name).CCR_zero_l_loading_std = std(all_data.(name).CCR_zero_l_loading(2:end));
    all_data.(name).CCR_zero_l_early_std = std(all_data.(name).CCR_zero_l_early(2:end));
    all_data.(name).CCR_zero_l_mid_std = std(all_data.(name).CCR_zero_l_mid(2:end));
    all_data.(name).CCR_zero_l_late_std = std(all_data.(name).CCR_zero_l_late(2:end));
    
    else
        none = 1;
    end
%%
end
