%% Combinging Gait Parametersparameters{j}
% Combine gait parameters across trials with the same condition
function [gait_metrics, gait_metrics_combined1] = gait_parameters_combined(gait_parameters)
%%
% Combine all variables besides speed
trials = fieldnames(gait_parameters);
sides = fieldnames(gait_parameters.(trials{1}));
parameters = fieldnames(gait_parameters.(trials{1}).(sides{1}));
for i = 1:length((sides))
    for j = 1:length(parameters)
        gait_metrics.control.(sides{i}).(parameters{j}) = [];
        gait_metrics.keeogo.(sides{i}).(parameters{j}) = [];

    end
end
%%
for i = 1:length(fieldnames(gait_parameters))
    for j = 1:length(sides)
        for k = 1:length(parameters)
            if contains(lower(trials(i)),'control')
                gait_metrics_combined = gait_parameters.(trials{i}).(sides{j}).(parameters{k});
                gait_metrics.control.(sides{j}).(parameters{k}) = [gait_metrics.control.(sides{j}).(parameters{k}), gait_metrics_combined(1:end)];
                
            elseif contains(lower(trials(i)),'keeogo')
                gait_metrics_combined = gait_parameters.(trials{i}).(sides{j}).(parameters{k});
                gait_metrics.keeogo.(sides{j}).(parameters{k}) = [gait_metrics.keeogo.(sides{j}).(parameters{k}), gait_metrics_combined(1:end)];
            end
        end
    end
end

for i = 1:length(sides)
    for j = 1:length(parameters)
        gait_metrics_combined1.control.(sides{i}).(parameters{j}).mean = mean(gait_metrics.control.(sides{i}).(parameters{j}),"omitnan");
        gait_metrics_combined1.control.(sides{i}).(parameters{j}).std = std(gait_metrics.control.(sides{i}).(parameters{j}),"omitnan");
        gait_metrics_combined1.keeogo.(sides{i}).(parameters{j}).mean = mean(gait_metrics.keeogo.(sides{i}).(parameters{j}),"omitnan");
        gait_metrics_combined1.keeogo.(sides{i}).(parameters{j}).std = std(gait_metrics.keeogo.(sides{i}).(parameters{j}),"omitnan");
    end
end

%% Speed 
gait_metrics.control.total.speed = [];
gait_metrics.keeogo.total.speed = [];
for i = 1:length(trials)
    if contains(lower(trials(i)),'control')
        gait_metrics_combined = gait_parameters.(trials{i}).total.speed;
        gait_metrics.control.total.speed = [gait_metrics.control.total.speed, gait_metrics_combined(1:end)];

    elseif contains(lower(trials(i)),'keeogo')
        gait_metrics_combined = gait_parameters.(trials{i}).total.speed;
        gait_metrics.keeogo.total.speed = [gait_metrics.keeogo.total.speed, gait_metrics_combined(1:end)];
    end
end

gait_metrics_combined1.control.total.speed.mean = mean(gait_metrics.control.total.speed,"omitnan");
gait_metrics_combined1.control.total.speed.std = std(gait_metrics.control.total.speed,"omitnan");
gait_metrics_combined1.keeogo.total.speed.mean = mean(gait_metrics.keeogo.total.speed,"omitnan");
gait_metrics_combined1.keeogo.total.speed.std = std(gait_metrics.keeogo.total.speed,"omitnan");


