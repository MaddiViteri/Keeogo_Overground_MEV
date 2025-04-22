%Herl-Keeogo Overground Protocol
%Author: Minori Iizuka
%Maintained/edited by: Maddi Viteri
%Last Updated: 04/17/2025 (MATLAB R2024a)(Opensim 4.5)

%-------------------------------------------------------------------------%

clear
import org.opensim.modeling.*
%% Copy data address

trial_path = uigetdir;%
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
disp(class(platform_files));

platform_files = platform_files(3:end);

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

% Mean w/o EMG
[gait_parameters_combined,loading_stance,early_stance,mid_stance,late_stance,full_stance,ROM,more_severe,less_severe] =...
    stats_table_no_emg_mean(gait_metrics_mean,...
    main_parameters_mean,participant_number,trial_name,severe,save_file_folder);
