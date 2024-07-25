function [gait_emg_metrics, gait_emg_metrics_combined1,gait_emg_metrics_trials]...
    = emg_parameters_combined_rms(all_data_norm2_rms)
% Does not include the swing and stance emg
%%
trials = fieldnames(all_data_norm2_rms);

% Get only emg trials that correspond with Vicon data
for i = 1:length(trials)
    aa = fieldnames(all_data_norm2_rms.(trials{i}));
    for j = 1:length(aa)
        if length(fieldnames(all_data_norm2_rms.(trials{i}))) <= 3
            disp('');
        else
            emg_gait_combined.(trials{i}).(aa{j}) = [];
        end
    end
end
%% Get actual trials and parameters that are used
used_trials = fieldnames(emg_gait_combined);

for i = 1:length(used_trials)
    if ~contains(lower(used_trials{i}), 'mvc') % skip fieldnames that contain 'mvc'
        used_parameters = fieldnames(emg_gait_combined.(used_trials{i}));
        % Now you can process 'substruct' as needed
        break; % Exit the loop after finding the first match
    end
end
% Initialize EMG gait meterics
for j = 1:length(used_parameters)
    gait_emg_metrics.control.(used_parameters{j}).values = [];
    gait_emg_metrics.keeogo.(used_parameters{j}).values = [];
end
%% Combine data

for i = 1:length(used_trials)
    % Don't analyze MVC trials
    if ~contains(lower(used_trials{i}),'mvc')
        % Only combine parameter data
        for j = 1:length(used_parameters)
            if ~isa(all_data_norm2_rms.(used_trials{i}).(used_parameters{j}),'struct')
                % Combine control trials
                if contains(lower(used_trials{i}),'control')
                    gait_emg_metrics_combined = all_data_norm2_rms.(used_trials{i}).(used_parameters{j});
                    gait_emg_metrics.control.(used_parameters{j}).values = [gait_emg_metrics.control.(used_parameters{j}).values, gait_emg_metrics_combined(2:end)];
                    gait_emg_metrics_trials.(used_trials{i}).(used_parameters{j})= gait_emg_metrics_combined;
                else % Combine Keeogo trials
                    gait_emg_metrics_combined = all_data_norm2_rms.(used_trials{i}).(used_parameters{j});
                    gait_emg_metrics.keeogo.(used_parameters{j}).values = [gait_emg_metrics.keeogo.(used_parameters{j}).values, gait_emg_metrics_combined(2:end)];
                    gait_emg_metrics_trials.(used_trials{i}).(used_parameters{j})= gait_emg_metrics_combined;
                end
            end
        end
    else
        disp('')
    end
end
%%
% Mean and StDev for control and keeogo
for j = 1:length(used_parameters)
        gait_emg_metrics_combined1.control.(used_parameters{j}).mean = mean(gait_emg_metrics.control.(used_parameters{j}).values,"omitnan");
        gait_emg_metrics_combined1.control.(used_parameters{j}).std = std(gait_emg_metrics.control.(used_parameters{j}).values,"omitnan");
        gait_emg_metrics_combined1.keeogo.(used_parameters{j}).mean = mean(gait_emg_metrics.keeogo.(used_parameters{j}).values,"omitnan");
        gait_emg_metrics_combined1.keeogo.(used_parameters{j}).std = std(gait_emg_metrics.keeogo.(used_parameters{j}).values,"omitnan");
end
end
