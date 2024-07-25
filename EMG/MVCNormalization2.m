% function [all_data] = MVCNormalization2(all_data,MVC,trial_names,muscle_names)
function [all_data] = MVCNormalization2(all_data,trial_names,actual_trials,muscle_names,type)

    %%
    for i = 1:length(trial_names)
        name = trial_names{i};
        for m = 1:length(muscle_names)
            muscle = muscle_names{m};
            if contains(lower(trial_names{i}),'control') && contains(lower(trial_names{i}),actual_trials)
                control_max_all_trials.(muscle)(i) = max(all_data.(name).(type).(muscle));
            end
        end
    end
    for  m = 1:length(muscle_names)
        muscle = muscle_names{m};
        control_max_all_trials.(muscle) = nonzeros(control_max_all_trials.(muscle));
        control_max.(muscle) = mean(control_max_all_trials.(muscle));
    end
    %%
    for t = 1:length(trial_names)
        name = trial_names{t};
        for m = 1:length(muscle_names)
            muscle = muscle_names{m};
            if isstruct(all_data.MVC)
                a = fieldnames(all_data.MVC);
                all_data.MVC_values.(muscle) = calculateMVC(all_data.MVC.(a{m}).(muscle));
                if max(all_data.(name).(type).(muscle)) < all_data.MVC_values.(muscle)
                    all_data.(name).(type).(muscle) = all_data.(name).(type).(muscle)./all_data.MVC_values.(muscle);
                else
                    all_data.(name).(type).(muscle) = all_data.(name).(type).(muscle)./control_max.(muscle);

                end
            else
                all_data.(name).(type).(muscle) = all_data.(name).(type).(muscle)./control_max.(muscle);
            end
        end
    end
    %%
end