function [emg_stance_metrics, emg_stance_metrics_combined, gait_emg_metrics_trials]...
    = emg_parameters_combined(all_data_norm2)
% Does not include the swing and stance emg
%%
trials = fieldnames(all_data_norm2);

% Get only emg trials that correspond with Vicon data
for i = 1:length(trials)
    aa = fieldnames(all_data_norm2.(trials{i}));
    for j = 1:length(aa)
        if length(fieldnames(all_data_norm2.(trials{i}))) <= 3
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
        all_used_parameters = fieldnames(emg_gait_combined.(used_trials{i}));
        for j = 1:numel(all_used_parameters)
            if length(all_data_norm2.(used_trials{i}).(all_used_parameters{j})) <= 1
                emg_gait_combined.(used_trials{i}) = ...
                    rmfield(emg_gait_combined.(used_trials{i}),(all_used_parameters{j}));
            end
        end
        used_parameters = fieldnames(emg_gait_combined.(used_trials{i}));
        break; % Exit the loop after finding the first match
    end
end
% Initialize EMG gait meterics
for j = 1:length(used_parameters)
    emg_stance_metrics.control.(used_parameters{j}) = [];
    emg_stance_metrics.keeogo.(used_parameters{j}) = [];
end
%% Combine data

for i = 1:length(used_trials)
    % Don't analyze MVC trials
    if ~contains(lower(used_trials{i}),'mvc')
        for j = 1:length(used_parameters)
            % Only combine parameter data
            if ~isa(all_data_norm2.(used_trials{i}).(used_parameters{j}),'struct') && length(all_data_norm2.(used_trials{i}).(used_parameters{j}))>1
                % Combine control trials
                if contains(lower(used_trials{i}),'control')
                    gait_emg_metrics_combined = all_data_norm2.(used_trials{i}).(used_parameters{j});
                    emg_stance_metrics.control.(used_parameters{j}) = [emg_stance_metrics.control.(used_parameters{j}), gait_emg_metrics_combined(2:end)];
                    gait_emg_metrics_trials.(used_trials{i}).(used_parameters{j})= all_data_norm2.(used_trials{i}).(used_parameters{j});
                else % Combine Keeogo trials
                    gait_emg_metrics_combined = all_data_norm2.(used_trials{i}).(used_parameters{j});
                    emg_stance_metrics.keeogo.(used_parameters{j}) = [emg_stance_metrics.keeogo.(used_parameters{j}), gait_emg_metrics_combined(2:end)];
                    gait_emg_metrics_trials.(used_trials{i}).(used_parameters{j})= all_data_norm2.(used_trials{i}).(used_parameters{j});
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
    emg_stance_metrics_combined.control.(used_parameters{j}).mean = mean(emg_stance_metrics.control.(used_parameters{j}),"omitnan");
    emg_stance_metrics_combined.control.(used_parameters{j}).std = std(emg_stance_metrics.control.(used_parameters{j}),"omitnan");
    emg_stance_metrics_combined.keeogo.(used_parameters{j}).mean = mean(emg_stance_metrics.keeogo.(used_parameters{j}),"omitnan");
    emg_stance_metrics_combined.keeogo.(used_parameters{j}).std = std(emg_stance_metrics.keeogo.(used_parameters{j}),"omitnan");
end
end
