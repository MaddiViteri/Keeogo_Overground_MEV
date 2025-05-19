% Herl-Keeogo Overground Protocol
% Author: Minori Iizuka
% Maintained/edited by: Maddi Viteri
% Last Updated: 05/02/2025 (MATLAB R2024a)(OpenSim 4.5)

% Ask user if they want to load new data or use workspace
choice = questdlg('Load new trial data or use existing workspace variables?', ...
    'Data Loading Option', ...
    'Load New', 'Use Existing', 'Load New');

use_existing = strcmp(choice, 'Use Existing');

% Clear only if user selected "Load New"
if ~use_existing
    clearvars -except use_existing
    clc;
end

import org.opensim.modeling.*

if ~use_existing
    %% Select Trial Folder
    trial_path = uigetdir;
    
    % ... all the data loading and processing code ...

    % Get list of all files in trial directory
    file_list = dir(trial_path);
    file_names = {file_list.name};
        
    %% Input Anthropometrics
    prompts = {'Weight (kg):', 'Height (m):', ...
               'Gender (Male = 1, Female = 0):', 'Keeogo (Yes = 1, No = 0):', ...
               'More Severe (Left = 1, Right = 0):', 'Participant # (XX):'};
    defaults = {'1','1','1','1','1','1'};
    anthropometrics = str2double(inputdlg(prompts, 'Anthropometrics', 1, defaults));
    
    [mass, height, gender, keeogo, severe, p_num] = deal(anthropometrics(1), ...
        anthropometrics(2), anthropometrics(3), anthropometrics(4), ...
        anthropometrics(5), string(anthropometrics(6)));
    
    %% Initialize Data Containers
    trial = struct(); 
    static = struct(); 
    knee = struct(); 
    ankle = struct();
    platform = struct(); 
    exo = struct();
    trial_name = {}; 
    static_trial_name = {};
    processed_files = {}; 
    no_force_data = {}; 
        
    opts = delimitedTextImportOptions("NumVariables", 7, ...
        "DataLines", [2, Inf], ...
        "Delimiter", "\t", ...
        "VariableNames", ["Time", "LeftBeltSpeed", "LeftBeltDistance", ...
                          "RightBeltSpeed", "RightBeltDistance", ...
                          "TreadmillPitch", "PlatformRotation"], ...
        "VariableTypes", repmat("double", 1, 7), ...
        "ExtraColumnsRule", "ignore", ...
        "EmptyLineRule", "read");
end
        
%% Parse Files
for n = 1:length(file_names)
    name = file_names{n};
    path = fullfile(trial_path, name);
    [~, base, ext] = fileparts(name);
    
    if strcmpi(ext, '.c3d')
        fprintf('Parsing file %d of %d: %s ...\n', n, length(file_names), name);
        processed_files{end+1} = name;

        try
            % Load and rotate C3D data
            c3d_data = osimC3D(path, 1);
            c3d_data.rotateData('x', -90);
            [markers, forces] = c3d_data.getAsStructs();

            % Check for missing force data
            if isempty(forces)
                warning('  No force data found in %s', name);
                no_force_data{end+1} = name;
            end

            % Categorize trial
            if contains(lower(base), 'static')
                static.(base).markers = markers;
                static.(base).force = forces;
                static_trial_name{end+1} = base;
            else
                trial.(base).markers = markers;
                trial.(base).force = forces;
                trial_name{end+1} = base;
            end

            % Construct path to treadmill-caren folder (sibling to 'Vicon/neutral')
            belt_folder = fullfile(fileparts(fileparts(trial_path)), 'treadmill-caren');
            belt_file = fullfile(belt_folder, [base, '.txt']);

            % Try to load .txt file if available
            if isfile(belt_file)
            try
                platform.(base) = readtable(belt_file, opts);
            catch readErr
                warning('  Failed to read %s: %s', belt_file, readErr.message);
            end
            else
                warning('  Belt data file not found for trial: %s', base);
            end

            fprintf('  Successfully parsed %s.\n', name);

        catch parseErr
            warning('  Failed to parse %s: %s', name, parseErr.message);
        end
    end
end


% Final checks
if isempty(static_trial_name) && isempty(trial_name)
    error('No C3D trials were found or processed.');
end

fprintf('Processed the following .c3d files:\n');
disp(processed_files');

if ~isempty(no_force_data)
    warning('The following files had no force data:');
    disp(no_force_data');
end

%% Process Kinetics
zcr_threshold = 50;
tic;

for i = 1:numel(trial_name)
    t_name = trial_name{i};
    markerStruct = trial.(t_name).markers;
    forceStruct = trial.(t_name).force;
    
    for j = 1:numel(static_trial_name)
        s_name = static_trial_name{j};
        static_markers = static.(s_name).markers;
        
        [K, A, markerStruct, forceStruct] = Inverse_3D_Dynamics_COM_7_18_debug( ...
            markerStruct, forceStruct, static_markers, anthropometrics, t_name);

        knee.(t_name) = K;
        ankle.(t_name) = A;
        trial.(t_name).markers = markerStruct;
        trial.(t_name).force = forceStruct;

        knee.(t_name) = Step_Cycle_MI_combined_brace_v_keeogo_24_03_13( ...
            forceStruct, K, height, mass, zcr_threshold);
    end
end
toc;



%% Gait Parameter Analysis
trial_fields = fieldnames(trial);
platform_fields = fieldnames(platform);
gait_parameters = struct();
main_parameters = struct();

disp('Available platform fields:');
disp(platform_fields);

for i = 1:numel(trial_fields)
    fname = trial_fields{i};
    
    disp(['Processing gait parameters for trial: ', fname]);

    if ~isfield(platform, fname)
        warning('No matching platform data for trial: %s. Skipping gait analysis.', fname);
        continue; % Skip to next trial
    end

    try
        gait_parameters.(fname) = SpatioTemporal( ...
            trial.(fname).force, ...
            trial.(fname).markers, ...
            platform.(fname), ...
            zcr_threshold);
    catch ME
        if contains(ME.message, 'Unrecognized field name "single_support')
            warning('SpatioTemporal failed for trial %s: %s. Skipping.', fname, ME.message);
            continue; % Skip to next trial
        else
            rethrow(ME); % Re-throw if it's a different error
        end
    end

    main_parameters.(fname) = stride_peak_values(trial.(fname).force, knee.(fname), zcr_threshold, fname);
end


% Check if any trials in gait_parameters have data
if ~isempty(fieldnames(gait_parameters))
    has_data = false;
    
    % Loop through trials and check for any valid data
    trial_fields = fieldnames(gait_parameters);
    for i = 1:numel(trial_fields)
        if isstruct(gait_parameters.(trial_fields{i})) && ~isempty(fieldnames(gait_parameters.(trial_fields{i})))
            has_data = true;
            break;
        end
    end

    if has_data
        %% Combine and Plot Results
        [gait_metrics_full, gait_metrics_mean] = gait_parameters_combined_debug(gait_parameters);
        [main_parameters_full, main_parameters_mean] = main_parameters_combined(main_parameters);

        Knee_combine = plot_function(knee);

        [gait_parameters_combined, loading_stance, early_stance, mid_stance, late_stance, full_stance, ROM, more_severe, less_severe] = stats_table_no_emg_mean(gait_metrics_mean, main_parameters_mean, p_num, trial_name, severe, 'stats_table_no_emg_mean');
    else
        error('No gait parameter data found. Please verify platform and trial data were processed.');
    end
end

disp('runtime complete');