%% test code
% emg_freq = 2000;
% 
% time = trial.neutral_control06.force.time;
% ForceLT = trial.neutral_control06.force.f1(:,2);
% ForceRT = trial.neutral_control06.force.f2(:,2);
% % 
% % ForceLT = sqrt(sum(ForceLT.^2,2));
% % ForceRT = sqrt(sum(ForceRT.^2,2)); 
% name = emg_trial_names{2};
% emg_time_full = all_data.(name).time;
% thresholds = 50;
% prep_phase = emg_freq*0.1;
function [all_data_norm2_rms] = swing_data_8_25_rms(all_data_norm2_rms,emg_trial_names,muscle_names,forceStruct,trial_name,thresholds)

emg_freq = 2000;
time = forceStruct.time;

ForceLT = forceStruct.f1(:,2);
ForceRT = forceStruct.f2(:,2);
prep_phase = emg_freq*0.1;
name = emg_trial_names;
emg_time_full = all_data_norm2_rms.(name).time;

%%

if contains(lower(name),lower(trial_name))
    m = 1:length(muscle_names);
    rquad = all_data_norm2_rms.(name).rms.(muscle_names{6});
    rham = all_data_norm2_rms.(name).rms.(muscle_names{4});
    rgas = all_data_norm2_rms.(name).rms.(muscle_names{2});
    lquad = all_data_norm2_rms.(name).rms.(muscle_names{5});
    lham = all_data_norm2_rms.(name).rms.(muscle_names{3});
    lgas = all_data_norm2_rms.(name).rms.(muscle_names{1});
    % muscle = muscle_names{m};            

    %% Right Leg
    [~, heel_strike_data,toe_off_data] = detect_cycles(ForceRT,time,thresholds);
    heel_time_r = cell2mat(heel_strike_data(:,1));
    toe_time_r = cell2mat(toe_off_data(:,1));
    nn_r = min(numel(heel_time_r),numel(toe_time_r));
    for time_index = 1:nn_r
        if heel_time_r(time_index) > toe_time_r(time_index)
            for ci = 2:nn_r
                emg_time = toe_time_r(ci):1/2000:heel_time_r(ci);
                emg_swing_range = round(round(toe_time_r(ci),3)*emg_freq)+60:round(round(heel_time_r(ci),3)*emg_freq)+60;
                for di = min(numel(emg_time),numel(emg_swing_range))
                    % Swing phase of each muscle
                    all_data_norm2_rms.(name).swing.rquad_swing(ci) = mean(rquad(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.rham_swing(ci) = mean(rham(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.rgas_swing(ci) = mean(rgas(emg_swing_range(1:di)))*100;
                    % Prep phase
                    all_data_norm2_rms.(name).swing.rquad_stride_prep(ci) = mean(rquad(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.rham_stride_prep(ci) = mean(rham(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.rgas_stride_prep(ci) = mean(rgas(emg_swing_range(round(di-prep_phase):di)))*100;
                end
                new_t1(1:100,ci) = linspace(toe_time_r(ci),heel_time_r(ci),100);
                all_data_norm2_rms.(name).rms_swing.rgas(:,ci) = interp1(emg_time_full,rgas,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.rham(:,ci) = interp1(emg_time_full,rham,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.rquad(:,ci) = interp1(emg_time_full,rquad,new_t1(:,ci)');
            end
        else
            for ci = 2:nn_r
                emg_time = toe_time_r(ci-1):1/2000:heel_time_r(ci);
                emg_swing_range = round(round(toe_time_r(ci-1),3)*emg_freq):round(round(heel_time_r(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_swing_range))
                    % Swing phase of each muscle
                    all_data_norm2_rms.(name).swing.rquad_swing(ci) = mean(rquad(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.rham_swing(ci) = mean(rham(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.rgas_swing(ci) = mean(rgas(emg_swing_range(1:di)))*100;
                    % Prep phase
                    all_data_norm2_rms.(name).swing.rquad_stride_prep(ci) = mean(rquad(emg_swing_range((di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.rham_stride_prep(ci) = mean(rham(emg_swing_range((di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.rgas_stride_prep(ci) = mean(rgas(emg_swing_range((di-prep_phase):di)))*100;
                end
                new_t1(1:100,ci) = linspace(toe_time_r(ci-1),heel_time_r(ci),100);
                all_data_norm2_rms.(name).rms_swing.rgas(:,ci) = interp1(emg_time_full,rgas,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.rham(:,ci) = interp1(emg_time_full,rham,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.rquad(:,ci) = interp1(emg_time_full,rquad,new_t1(:,ci)');
            end
        end
    end

    %% Left Leg
    [~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds);
    heel_time_l = cell2mat(heel_strike_data_l(:,1));
    toe_time_l = cell2mat(toe_off_data_l(:,1));
    nn_l = min(numel(heel_time_l),numel(toe_time_l));
    for time_index = 1:nn_l
        if heel_time_l(time_index) > toe_time_l(time_index)
            for ci = 1:nn_l
                emg_time = toe_time_l(ci):1/2000:heel_time_l(ci);
                emg_swing_range = round(round(toe_time_l(ci),3)*emg_freq):round(round(heel_time_l(ci),3)*emg_freq);
                for di = min(numel(emg_time),numel(emg_swing_range))
                    % Swing phase of each muscle
                    all_data_norm2_rms.(name).swing.lquad_swing(ci) = mean(lquad(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.lham_swing(ci) = mean(lham(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.lgas_swing(ci) = mean(lgas(emg_swing_range(1:di)))*100;
                    % Prep phase
                    all_data_norm2_rms.(name).swing.lquad_stride_prep(ci) = mean(lquad(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.lham_stride_prep(ci) = mean(lham(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.lgas_stride_prep(ci) = mean(lgas(emg_swing_range(round(di-prep_phase):di)))*100;
                end
                new_t1(1:100,ci) = linspace(toe_time_l(ci),heel_time_l(ci),100);
                all_data_norm2_rms.(name).rms_swing.lgas(:,ci) = interp1(emg_time_full,lgas,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.lham(:,ci) = interp1(emg_time_full,lham,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.lquad(:,ci) = interp1(emg_time_full,lquad,new_t1(:,ci)');
            end
        else
            for ci = 2:nn_l
                emg_time = toe_time_l(ci-1):1/2000:heel_time_l(ci);
                emg_swing_range = round(round(toe_time_l(ci-1),3)*emg_freq)+60:round(round(heel_time_l(ci),3)*emg_freq)+60;
                for di = min(numel(emg_time),numel(emg_swing_range))
                    % Swing phase of each muscle
                    all_data_norm2_rms.(name).swing.lquad_swing(ci) = mean(lquad(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.lham_swing(ci) = mean(lham(emg_swing_range(1:di)))*100;
                    all_data_norm2_rms.(name).swing.lgas_swing(ci) = mean(lgas(emg_swing_range(1:di)))*100;
                    % Prep phase
                    all_data_norm2_rms.(name).swing.lquad_stride_prep(ci) = mean(lquad(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.lham_stride_prep(ci) = mean(lham(emg_swing_range(round(di-prep_phase):di)))*100;
                    all_data_norm2_rms.(name).swing.lgas_stride_prep(ci) = mean(lgas(emg_swing_range(round(di-prep_phase):di)))*100;
                end
                new_t1(1:100,ci) = linspace(toe_time_l(ci-1),heel_time_l(ci),100);
                all_data_norm2_rms.(name).rms_swing.lgas(:,ci) = interp1(emg_time_full,lgas,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.lham(:,ci) = interp1(emg_time_full,lham,new_t1(:,ci)');
                all_data_norm2_rms.(name).rms_swing.lquad(:,ci) = interp1(emg_time_full,lquad,new_t1(:,ci)');
            end
        end
    end

else
    none =0;
end
%%
end
