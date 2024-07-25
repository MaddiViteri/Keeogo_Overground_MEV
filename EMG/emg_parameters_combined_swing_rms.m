function [gait_emg_swing_metrics, gait_emg_swing_metrics_combined,used_trials] ...
    = emg_parameters_combined_swing_rms(all_data_norm2_rms)
% f
%%
trials = fieldnames(all_data_norm2_rms);

for i = 1:length(trials)
    aa = fieldnames(all_data_norm2_rms.(trials{i}));
    for j = 1:length(aa)
        if length(fieldnames(all_data_norm2_rms.(trials{i}))) <= 3
            disp('');
        else
            emg_gait_combined.(trials{i}).(aa{j}).swing = [];
        end
    end
end
%% Get actual trials and parameters that are used
used_trials = fieldnames(emg_gait_combined);

for i = 1:length(used_trials)
    if ~contains(lower(used_trials{i}), 'mvc') % skip fieldnames that contain 'mvc'
        all_used_parameters = fieldnames(all_data_norm2_rms.(used_trials{i}).swing);
        for j = 1:numel(all_used_parameters)
            if length(all_data_norm2_rms.(used_trials{i}).swing.(all_used_parameters{j})) <= 1
                emg_gait_combined.(used_trials{i}).swing = ...
                    rmfield(emg_gait_combined.(used_trials{i}).swing,(all_used_parameters{j}));
            end
        end
        used_parameters = fieldnames(all_data_norm2_rms.(used_trials{i}).swing);
        break; % Exit the loop after finding the first match
    end
end

for j = 1:length(used_parameters)
    gait_emg_swing_metrics.control.(used_parameters{j}) = [];
    gait_emg_swing_metrics.keeogo.(used_parameters{j}) = [];
end
%% Combine data

for i = 1:length(used_trials)
    % Don't analyze MVC trials
    if ~contains(lower(used_trials{i}),'mvc')
        for j = 1:length(used_parameters)
            % Combine control trials
            if ~isa(all_data_norm2_rms.(used_trials{i}).swing.(used_parameters{j}),'struct')
                if contains(lower(used_trials{i}),'control')
                    gait_emg_metrics_trial = all_data_norm2_rms.(used_trials{i}).swing.(used_parameters{j});
                    gait_emg_swing_metrics.control.(used_parameters{j}) = [gait_emg_swing_metrics.control.(used_parameters{j}), gait_emg_metrics_trial(2:end)];
                    
                else % Combine Keeogo trials
                    gait_emg_metrics_trial = all_data_norm2_rms.(used_trials{i}).swing.(used_parameters{j});
                    gait_emg_swing_metrics.keeogo.(used_parameters{j}) = [gait_emg_swing_metrics.keeogo.(used_parameters{j}), gait_emg_metrics_trial(2:end)];
                end
            end
        end
    else
        disp('')
    end
end
%%
% Mean and StDev
for j = 1:length(used_parameters)
    gait_emg_swing_metrics_combined.control.(used_parameters{j}).mean = mean(gait_emg_swing_metrics.control.(used_parameters{j}),"omitnan");
    gait_emg_swing_metrics_combined.control.(used_parameters{j}).std = std(gait_emg_swing_metrics.control.(used_parameters{j}),"omitnan");
    gait_emg_swing_metrics_combined.keeogo.(used_parameters{j}).mean = mean(gait_emg_swing_metrics.keeogo.(used_parameters{j}),"omitnan");
    gait_emg_swing_metrics_combined.keeogo.(used_parameters{j}).std = std(gait_emg_swing_metrics.keeogo.(used_parameters{j}),"omitnan");
end

end
