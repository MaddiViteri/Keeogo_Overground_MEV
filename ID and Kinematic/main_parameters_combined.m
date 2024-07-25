%% Combinging Gait Parameters
% Combine gait parameters across trials with the same condition

function [main_metrics, main_parameters_combined1] = main_parameters_combined(main_parameters)
trials = fieldnames(main_parameters);
sides = fieldnames(main_parameters.(trials{1}));
parameters = fieldnames(main_parameters.(trials{1}).(sides{1}));

% Initialize structure
for x = 1:length(trials)
    for y = 1:length(sides)
        for z = 1:length(parameters)
            if ~isa(main_parameters.(trials{x}).(sides{y}).(parameters{z}),'struct')
                main_metrics.control.(sides{y}).(parameters{z}) = [];
                main_metrics.keeogo.(sides{y}).(parameters{z}) = [];
            else
                phase = fieldnames(main_parameters.(trials{x}).(sides{y}).(parameters{z}));
                for k = 1:length(phase)
                    main_metrics.control.(sides{y}).(parameters{z}).(phase{k}) = [];
                    main_metrics.keeogo.(sides{y}).(parameters{z}).(phase{k}) = [];
                end
            end
        end
    end
end

% Assign variables
for x = 1:length(trials)
    for y = 1:length(sides)
        for z = 1:length(parameters)
            if ~isa(main_parameters.(trials{x}).(sides{y}).(parameters{z}),'struct')
                if contains(lower(trials(x)),'control')
                    main_parameters_combined = main_parameters.(trials{x}).(sides{y}).(parameters{z});
                    main_metrics.control.(sides{y}).(parameters{z}) =...
                        [main_metrics.control.(sides{y}).(parameters{z}),...
                        main_parameters_combined];
                    
                elseif contains(lower(trials(x)),'keeogo')
                    main_parameters_combined = main_parameters.(trials{x}).(sides{y}).(parameters{z});
                    main_metrics.keeogo.(sides{y}).(parameters{z}) =...
                        [main_metrics.keeogo.(sides{y}).(parameters{z}),...
                        main_parameters_combined];
                end
            else
                for k = 1:length(phase)
                    if contains(lower(trials(x)),'control')
                        main_parameters_combined = main_parameters.(trials{x}).(sides{y}).(parameters{z}).(phase{k});
                        main_metrics.control.(sides{y}).(parameters{z}).(phase{k}) = ...
                            [main_metrics.control.(sides{y}).(parameters{z}).(phase{k}),...
                            main_parameters_combined];
                    
                    elseif contains(lower(trials(x)),'keeogo')
                        main_parameters_combined = main_parameters.(trials{x}).(sides{y}).(parameters{z}).(phase{k});
                        main_metrics.keeogo.(sides{y}).(parameters{z}).(phase{k}) =...
                            [main_metrics.keeogo.(sides{y}).(parameters{z}).(phase{k}),...
                            main_parameters_combined];
                    end
                end
            end
        end
    end
end

% Mean and StDev Calculated for Parameters
for x = 1:length(trials)
    for y = 1:length(sides)
        for z = 1:length(parameters)
            if ~isa(main_parameters.(trials{x}).(sides{y}).(parameters{z}),'struct')
                main_parameters_combined1.control.(sides{y}).(parameters{z}).mean = mean(main_metrics.control.(sides{y}).(parameters{z}),"omitnan");
                main_parameters_combined1.control.(sides{y}).(parameters{z}).std = std(main_metrics.control.(sides{y}).(parameters{z}),"omitnan");
                main_parameters_combined1.keeogo.(sides{y}).(parameters{z}).mean = mean(main_metrics.keeogo.(sides{y}).(parameters{z}),"omitnan");
                main_parameters_combined1.keeogo.(sides{y}).(parameters{z}).std = std(main_metrics.keeogo.(sides{y}).(parameters{z}),"omitnan");
            else
                for k = 1:length(phase)
                    main_parameters_combined1.control.(sides{y}).(parameters{z}).(phase{k}).mean = ...
                    mean(main_metrics.control.(sides{y}).(parameters{z}).(phase{k}),"omitnan");
                    main_parameters_combined1.control.(sides{y}).(parameters{z}).(phase{k}).std = std(main_metrics.control.(sides{y}).(parameters{z}).(phase{k}),"omitnan");
                    main_parameters_combined1.keeogo.(sides{y}).(parameters{z}).(phase{k}).mean = mean(main_metrics.keeogo.(sides{y}).(parameters{z}).(phase{k}),"omitnan");
                    main_parameters_combined1.keeogo.(sides{y}).(parameters{z}).(phase{k}).std = std(main_metrics.keeogo.(sides{y}).(parameters{z}).(phase{k}),"omitnan");
                end
            end
        end
    end
end
