% Debug script to identify marker issues
clear;
import org.opensim.modeling.*

% Get the trial path
trial_path = uigetdir;

% Find static trial files
file_names = {dir(trial_path).name};
static_file_path = '';

for n = 3:length(file_names)
    file_name = file_names{n};
    full_path = fullfile(trial_path, file_name);
    
    % Check if it's a static C3D file
    if endsWith(lower(file_name), '.c3d') && contains(lower(file_name), 'static')
        static_file_path = full_path;
        disp(['Found static file: ' file_name]);
        break;
    end
end

if isempty(static_file_path)
    error('No static file found');
end

% Load the static file
c3d = osimC3D(static_file_path, 1);
c3d.rotateData('x', -90);
[markers, ~] = c3d.getAsStructs();

% Check for the problematic marker
if isfield(markers, 'LLatMal')
    disp('LLatMal marker is present in the data');
else
    disp('LLatMal marker is NOT present in the data');
    disp('Available markers:');
    disp(fieldnames(markers));
end

% Save the marker data for further inspection
save('debug_markers.mat', 'markers');

% Try a simpler direct test of the calibration function
try
    disp('Testing calibrationmatrix_6_27 with the loaded markers...');
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = calibrationmatrix_6_27(markers, markers);
    disp('Function executed successfully');
catch ME
    disp('Error occurred:');
    disp(ME.message);
    if contains(ME.message, 'LLatMal')
    end
end