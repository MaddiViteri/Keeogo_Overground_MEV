function [gait_metrics, gait_metrics_combined1] = gait_parameters_combined_debug(gait_parameters)
    % Initialize gait_metrics
    gait_metrics = struct();
    gait_metrics_combined1 = struct();
    
    % Check if gait_parameters is empty or doesn't have the expected structure
    if isempty(gait_parameters)
        error('gait_parameters is empty. Please check the input data.');
    end
    
    % Display the structure of gait_parameters for debugging purposes
    disp('Structure of gait_parameters:');
    disp(gait_parameters);
    
    % Extract trials and check if it contains trials
    trials = fieldnames(gait_parameters);
    
    if isempty(trials)
        error('No trials found in gait_parameters. Please ensure data is properly loaded.');
    end
    
    % Check if there is at least one side (e.g., 'left', 'right')
    sides = fieldnames(gait_parameters.(trials{1}));
    
    if isempty(sides)
        error('No sides found in gait_parameters for trial: %s.', trials{1});
    end
    
    % Check if there are parameters available (e.g., 'stride', 'step')
    parameters = fieldnames(gait_parameters.(trials{1}).(sides{1}));
    
    if isempty(parameters)
        error('No parameters found in gait_parameters for side: %s, trial: %s.', sides{1}, trials{1});
    end
    
    % Initialize gait_metrics for control and keeogo
    for i = 1:length(sides)
        for j = 1:length(parameters)
            gait_metrics.control.(sides{i}).(parameters{j}) = [];
            gait_metrics.keeogo.(sides{i}).(parameters{j}) = [];
        end
    end

    % Combine metrics for control and keeogo trials
    for i = 1:length(trials)
        for j = 1:length(sides)
            for k = 1:length(parameters)
                if contains(lower(trials{i}), 'control')
                    gait_metrics_combined = gait_parameters.(trials{i}).(sides{j}).(parameters{k});
                    gait_metrics.control.(sides{j}).(parameters{k}) = [gait_metrics.control.(sides{j}).(parameters{k}), gait_metrics_combined(1:end)];
                elseif contains(lower(trials{i}), 'keeogo')
                    gait_metrics_combined = gait_parameters.(trials{i}).(sides{j}).(parameters{k});
                    gait_metrics.keeogo.(sides{j}).(parameters{k}) = [gait_metrics.keeogo.(sides{j}).(parameters{k}), gait_metrics_combined(1:end)];
                end
            end
        end
    end

    % Calculate mean and std for control and keeogo metrics
    for i = 1:length(sides)
        for j = 1:length(parameters)
            gait_metrics_combined1.control.(sides{i}).(parameters{j}).mean = mean(gait_metrics.control.(sides{i}).(parameters{j}), 'omitnan');
            gait_metrics_combined1.control.(sides{i}).(parameters{j}).std = std(gait_metrics.control.(sides{i}).(parameters{j}), 'omitnan');
            gait_metrics_combined1.keeogo.(sides{i}).(parameters{j}).mean = mean(gait_metrics.keeogo.(sides{i}).(parameters{j}), 'omitnan');
            gait_metrics_combined1.keeogo.(sides{i}).(parameters{j}).std = std(gait_metrics.keeogo.(sides{i}).(parameters{j}), 'omitnan');
        end
    end

    %% Speed
    gait_metrics.control.total.speed = [];
    gait_metrics.keeogo.total.speed = [];
    for i = 1:length(trials)
        if contains(lower(trials{i}), 'control')
            gait_metrics_combined = gait_parameters.(trials{i}).total.speed;
            gait_metrics.control.total.speed = [gait_metrics.control.total.speed, gait_metrics_combined(1:end)];
        elseif contains(lower(trials{i}), 'keeogo')
            gait_metrics_combined = gait_parameters.(trials{i}).total.speed;
            gait_metrics.keeogo.total.speed = [gait_metrics.keeogo.total.speed, gait_metrics_combined(1:end)];
        end
    end

    gait_metrics_combined1.control.total.speed.mean = mean(gait_metrics.control.total.speed, 'omitnan');
    gait_metrics_combined1.control.total.speed.std = std(gait_metrics.control.total.speed, 'omitnan');
    gait_metrics_combined1.keeogo.total.speed.mean = mean(gait_metrics.keeogo.total.speed, 'omitnan');
    gait_metrics_combined1.keeogo.total.speed.std = std(gait_metrics.keeogo.total.speed, 'omitnan');
end
