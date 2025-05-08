clear
import org.opensim.modeling.*

% Use the same trial path as in your main script
trial_path = uigetdir;
disp(['Selected path: ' trial_path]);

% Get file names
file_names = {dir(trial_path).name};
disp('Files in directory:');
disp(file_names(3:end)); % Skip . and .. entries

% Find and process static trial files
for n = 3:length(file_names)
    file_name = file_names{n};
    file_path = fullfile(trial_path, file_name);
    
    % Check if it's a C3D file and contains 'static'
    if endsWith(lower(file_name), '.c3d') && contains(lower(file_name), 'static')
        disp(['Processing static file: ' file_name]);
        
        % Load the static trial
        static_trial = osimC3D(file_path, 1);
        static_trial.rotateData('x', -90);
        
        % Get marker data
        [markers, ~] = static_trial.getAsStructs();
        
        % Display all marker names
        disp('Available marker names:');
        all_markers = fieldnames(markers);
        disp(all_markers);
        
        % Save marker names to a file for reference
        save('marker_names.mat', 'all_markers');
        disp('Marker names saved to marker_names.mat');
        
        % Break after finding first static trial
        break;
    end
end

if ~exist('markers', 'var')
    disp('No static trials found in the selected directory.');
end