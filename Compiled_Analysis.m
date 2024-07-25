clear
import org.opensim.modeling.*
%% Copy data address

trial_path ='E:\New folder (2)\06\Downhill';%
% trial_path = 'F:\Keeogo_Exoskeleton\02\Downhill';


file_names = {dir(trial_path).name}; % Write each trial and static

% Input dialog for anthropometrics
anthropometrics = str2double(inputdlg({'Weight (kg):','Height (m):',...
    'Gender (Male = 1, Female = 0) :','Keeogo (Yes = 1, No = 0)',...
    'More Severe (Left = 1, Right = 0):','Participant # (XX):'},'Anthropometrics',1,{'1','1','1','1','1','1'}));
% Assign input anthropometrics to mass, height, gender, and participant #
mass = anthropometrics(1);
height = anthropometrics(2);
gender = anthropometrics(3);
keeogo = anthropometrics(4);
severe = anthropometrics(5);
participant_number = string(anthropometrics(6));
%% Get files
for n = 3:length(file_names) % First 2 files are blank in file names
    true_name{n} = file_names{n}; % reassign file names to true names
    trial_file_name = true_name{n}(1:end-4); % get just file name, not file type
    file_path = fullfile(trial_path,file_names); % append file name with to folder where data is located

    % Categorize data into trial, static trial, platform data, and Keeogo data
    if strcmpi(file_path{n}(end-3:end),'.c3d') % Trial and Static Data
        c3d_trial = file_path{n};
        
        if ~contains(lower(c3d_trial),'static') % Trial data
           trial_name{n} = true_name{n}(1:end-4);
           c3d = osimC3D(c3d_trial,1); % Extract data
           c3d.rotateData('x',-90) % Rotate data

           % Assign structures for marker and force data
           [trial.(trial_name{n}).markers, trial.(trial_name{n}).force] =  c3d.getAsStructs();
           
        else % Static data
           static_trial_name{n} = true_name{n}(1:end-4);
           static_trial = osimC3D(c3d_trial,1); % Extract data
           static_trial.rotateData('x',-90) % Rotate data
           % Assign marker and force structure
           [static.(static_trial_name{n}).markers, static.(static_trial_name{n}).force] = static_trial.getAsStructs();
        end

    % Treadmill/Platform Data
    elseif strcmpi(file_path{n}(end-3:end),'.txt') 
        platform_files = file_path{n};
        platform_name{n} = true_name{n}(1:end-4);

        opts = delimitedTextImportOptions("NumVariables", 7);
        % Specify range and delimiter
        opts.DataLines = [2, Inf];
        opts.Delimiter = "\t";
        % Specify column names and types
        opts.VariableNames = ["Time", "LeftBeltSpeed", "LeftBeltDistance", "RightBeltSpeed", "RightBeltDistance", "TreadmillPitch", "PlatformRotation"];
        opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        % Data for Platforms
        gait_parameter_trial = platform_files;
        platform.(platform_name{n}) = table2array(readtable(strcat(file_path{n}), opts));
        
    % Exoskeleton Files
    elseif strcmpi(file_path{n}(end-3:end),'.txt') && contains(file_path{n}(1:end-4),'exo') 
        exo_files = file_path{n};
        exo_name{n} = true_name{n}(1:end-4);
        opts = delimitedTextImportOptions("NumVariables", 7);
        % Specify range and delimiter
        opts.DataLines = [2, Inf];
        opts.Delimiter = "\t";
        % Specify column names and types
        opts.VariableNames = ["Time", "LeftBeltSpeed", "LeftBeltDistance", "RightBeltSpeed", "RightBeltDistance", "TreadmillPitch", "PlatformRotation"];
        opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        exo_parameter_trial = exo_files;
        exo.(exo_name{n}) = table2array(readtable(strcat(file_path{n}), opts));
    end

end
%% Process Kinetic data
% Remove empty cells
trial_name = trial_name(~cellfun('isempty',trial_name)); 
static_trial_name = static_trial_name(~cellfun('isempty',static_trial_name));
platform_files = platform_files(1,3:end);
zcr_threshold = 50;
tic
% Process all trials for the Slope
for x = 1:length(trial_name)
    for y = 1:length(static_trial_name)
        % Run Inverse Dynamic Calculation Function 
        [knee.(trial_name{x}), ankle.(trial_name{x}),trial.(trial_name{x}).markers, trial.(trial_name{x}).force] = Inverse_3D_Dynamics_COM_7_18(trial.(trial_name{x}).markers, trial.(trial_name{x}).force, static.(static_trial_name{y}).markers,anthropometrics,trial_name{x});
        % [rknee.(trial_name{x}), lknee.(trial_name{x})] = Inverse_3D_Dynamics_GRF(trial.(trial_name{x}).markers, trial.(trial_name{x}).force, static.(static_trial_name{y}).markers,anthropometrics,trial_name{x});
        
        % Analyze kinetic variables for each step
        [knee.(trial_name{x})] = Step_Cycle_MI_combined_brace_v_keeogo_24_03_13(trial.(trial_name{x}).force,knee.(trial_name{x}),height,mass,zcr_threshold);
    end
end
toc
%% Process Gait Parameter Data
a = fieldnames(platform);

% Analyze Spatio-Temporal values
for i = 1:length(trial_name)
    [gait_parameters.(trial_name{i})] = SpatioTemporal(trial.(trial_name{i}).force,...
        trial.(trial_name{i}).markers,platform.(a{i}),zcr_threshold);
end

% Analyze Kinetic Variables and Main Parameters
for i = 1:length(trial_name)
     [main_parameters.(trial_name{i})] = stride_peak_values(trial.(trial_name{i}).force,...
         knee.(trial_name{i}),zcr_threshold,trial_name{i});
end

%% Combine data
% Combine the Keeogo and Control trial data separately for comparison
[gait_metrics_full, gait_metrics_mean] = gait_parameters_combined(gait_parameters);
[main_parameters_full, main_parameters_mean] = main_parameters_combined(main_parameters);
%% Plot Data of the Combined Trials
[Knee_combine] = plot_function(knee);

%% EMG
% emg_freq = 2000; 
% % Set Path for EMG data Location
% % emg_path = 'C:\Users\MII36\University of Pittsburgh\Keeogo Project - Documents\Pilot Testing\09-08-23\09-08-23\Keeogo_Exoskeleton\m, l\Comparison_setup';
% % emg_path = 'C:\Users\MII36\University of Pittsburgh\Keeogo Project - Documents\Pilot Testing\09-08-23\09-08-23\Keeogo_Exoskeleton\m, l\Our_set_up';
% emg_path = 'C:\Users\MII36\University of Pittsburgh\Keeogo Project - Documents\Data Processing\Pilot Testing\7-23-23\EMG\Keeogo_Exoskeleton\G, D';
% % emg_path = 'E:\New folder (2)\06\Keeogo_Exoskeleton\06';
% 
% % Assign muscle names
% muscle_names = {'lgas','rgas','lham','rham','lquad','rquad'};
% 
% trial_name = trial_name(~cellfun('isempty',trial_name));
% 
% % Filter EMG Data
% [all_data,emg_trial_names] = emgAnalysis(emg_path,muscle_names);
% 
% % Normalize EMG data for EMG and RMS
% [all_data_norm2] = MVCNormalization2(all_data,emg_trial_names,trial_name,muscle_names,'emg');
% [all_data_norm2_rms] = MVCNormalization2(all_data,emg_trial_names,trial_name,muscle_names,'rms');
% 
% % Process all trials for the Slope
% for x = 1:length(trial_name)
%     for y = 1:length(emg_trial_names)
% 
%         % Calculate for Co-activation (CCI) (directionality)
%         [all_data_norm2] = coactivationratio_zero(all_data_norm2,emg_trial_names{y},muscle_names,trial.(trial_name{x}).force,trial_name{x},zcr_threshold);
%         % Calculate Co-activation (CCR) (Magnitude)
%         [all_data_norm2_rms] = coactivationratio_09_01_rms(all_data_norm2_rms,emg_trial_names{y},muscle_names,trial.(trial_name{x}).force,trial_name{x},zcr_threshold);
% 
%         % Calculate swing values
%         [all_data_norm2_rms] = swing_data_8_25_rms(all_data_norm2_rms,emg_trial_names{y},muscle_names,trial.(trial_name{x}).force,trial_name{x},zcr_threshold);
%     end
% end

%% Combine EMG trials
% [emg_stance_metrics_rms, emg_stance_metrics_combined_rms,emg_metrics_trials_rms]...
%     = emg_parameters_combined(all_data_norm2_rms);
% 
% [emg_stance_metrics, emg_stance_metrics_combined,emg_metrics_trials]...
%     = emg_parameters_combined(all_data_norm2);
% 
% [emg_swing_metrics, emg_swing_metrics_combined,used_trials]...
%     = emg_parameters_combined_swing_rms(all_data_norm2_rms);

%% Plot EMG
% [stuff] = emg_plots(emg_stance_metrics_combined_rms,emg_stance_metrics_combined,emg_swing_metrics_combined);
%% Table
save_file_folder =  'C:\Users\MII36\University of Pittsburgh\Keeogo Project - Documents\Data Processing\Keeogo Data Analysis\Excel_Files\Downhill';

% % ID and EMG
% [stride_length,loading_stance,early_stance,mid_stance,late_stance,ROM] =...
%     stats_table(gait_metrics_full,gait_emg_metrics,gait_emg_metrics_rms,...
%     main_parameters_full,participant_number);
% [loading_stance,early_stance,mid_stance,late_stance,combined] =...
%     stats_table(gait_emg_metrics,gait_emg_metrics_rms,participant_number);
% [loading_stance,early_stance,mid_stance,late_stance] =...
%     stats_table(gait_emg_metrics,gait_emg_metrics_rms,participant_number);
% % No EMG
% [stride_length,loading_stance,early_stance,mid_stance,late_stance,full_stance,ROM] =...
%     stats_table_no_emg(gait_metrics_full,...
%     main_parameters_full,participant_number);
% Mean w/o EMG
[gait_parameters_combined,loading_stance,early_stance,mid_stance,late_stance,full_stance,ROM,more_severe,less_severe] =...
    stats_table_no_emg_mean(gait_metrics_mean,...
    main_parameters_mean,participant_number,trial_name,severe,save_file_folder);
% % Mean w/ EMG
% [gait_parameters_combined,loading_stance,early_stance,mid_stance,late_stance,full_stance,ROM,swing] =...
%     stats_table_w_emg_mean(gait_metrics_mean,...
%     main_parameters_mean,participant_number);
%% Stats analysis
% [loading_stance,mid_stance,late_stance,ROM,stride_length]=...
%     statistcal_tests(loading_stance,mid_stance,late_stance,...
%     gait_emg_metrics,gait_emg_metrics_rms,main_parameters_full,gait_metrics_full);
%   % [loading_stance,early_stance,mid_stance,late_stance]=...
%   %   statistcal_tests(loading_stance,early_stance,mid_stance,late_stance,gait_emg_metrics,gait_emg_metrics_rms);