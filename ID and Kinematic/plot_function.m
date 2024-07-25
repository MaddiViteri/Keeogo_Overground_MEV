function [Knee_combine] = plot_function(Knee)
% Initialize L/Knee.Righte_combine data
trials = fieldnames(Knee);
Knee_combine.Left.control.flexion_extension.values = [];
Knee_combine.Right.control.flexion_extension.values = [];
Knee_combine.Total.control.flexion_extension.values = [];
Knee_combine.Left.control.ab_adduction.values = [];
Knee_combine.Right.control.ab_adduction.values = [];
Knee_combine.Total.control.ab_adduction.values = [];
Knee_combine.Left.keeogo.flexion_extension.values= [];
Knee_combine.Right.keeogo.flexion_extension.values = [];
Knee_combine.Total.keeogo.flexion_extension.values = [];
Knee_combine.Left.keeogo.ab_adduction.values= [];
Knee_combine.Right.keeogo.ab_adduction.values = [];
Knee_combine.Total.keeogo.ab_adduction.values = [];

% Left Trials
for i = 1:length(trials)
    % Split between control and keeogo
    if contains(lower(trials(i)),'control')
        flex_ext_control_l = Knee.(trials{i}).Left.new_flexionmoment(:,2:end);
        ab_ad_duction_control_l = Knee.(trials{i}).Left.new_ab_add_moment(:,2:end);
        Knee_combine.Left.control.flexion_extension.values = [Knee_combine.Left.control.flexion_extension.values, flex_ext_control_l];
        Knee_combine.Left.control.ab_adduction.values = [Knee_combine.Left.control.ab_adduction.values, ab_ad_duction_control_l];
    elseif contains(lower(trials(i)),'keeogo')
        flex_ext_keeogo_l =  Knee.(trials{i}).Left.new_flexionmoment(:,2:end);
        ab_ad_duction_keeogo_l = Knee.(trials{i}).Left.new_ab_add_moment(:,2:end);
        Knee_combine.Left.keeogo.flexion_extension.values = [Knee_combine.Left.keeogo.flexion_extension.values, flex_ext_keeogo_l];
        Knee_combine.Left.keeogo.ab_adduction.values = [Knee_combine.Left.keeogo.ab_adduction.values, ab_ad_duction_keeogo_l];
    end
end

% Right Trials
for i = 1:length(trials)
    % Split between control and keeogo
    if contains(lower(trials(i)),'control')
        flex_ext_control_r = Knee.(trials{i}).Right.new_flexionmoment(:,2:end);
        ab_ad_duction_control_r = Knee.(trials{i}).Right.new_ab_add_moment(:,2:end);
        Knee_combine.Right.control.flexion_extension.values = [Knee_combine.Right.control.flexion_extension.values, flex_ext_control_r];
        Knee_combine.Right.control.ab_adduction.values = [Knee_combine.Right.control.ab_adduction.values, ab_ad_duction_control_r];
    elseif contains(lower(trials(i)),'keeogo')
        flex_ext_keeogo_r =  Knee.(trials{i}).Right.new_flexionmoment(:,2:end);
        ab_ad_duction_keeogo_r = Knee.(trials{i}).Right.new_ab_add_moment(:,2:end);
        Knee_combine.Right.keeogo.flexion_extension.values = [Knee_combine.Right.keeogo.flexion_extension.values, flex_ext_keeogo_r];
        Knee_combine.Right.keeogo.ab_adduction.values = [Knee_combine.Right.keeogo.ab_adduction.values, ab_ad_duction_keeogo_r];
    end
end

%Combined Right and Left Trials
for i = 1:length(trials)
    % Split between control and keeogo
    if contains(lower(trials(i)),'control')
        flex_ext_control_tot = [Knee.(trials{i}).Left.new_flexionmoment(:,2:end),...
            Knee.(trials{i}).Right.new_flexionmoment(:,2:end)];
        ab_ad_duction_control_tot = [-Knee.(trials{i}).Left.new_ab_add_moment(:,2:end),...
            Knee.(trials{i}).Right.new_ab_add_moment(:,2:end)];
        Knee_combine.Total.control.flexion_extension.values = [Knee_combine.Total.control.flexion_extension.values, flex_ext_control_tot];
        Knee_combine.Total.control.ab_adduction.values = [Knee_combine.Total.control.ab_adduction.values, ab_ad_duction_control_tot];
    elseif contains(lower(trials(i)),'keeogo')
        flex_ext_keeogo_tot =  [Knee.(trials{i}).Left.new_flexionmoment(:,2:end),...
            Knee.(trials{i}).Right.new_flexionmoment(:,2:end)];
        ab_ad_duction_keeogo_tot = [-Knee.(trials{i}).Left.new_ab_add_moment(:,2:end),...
            Knee.(trials{i}).Right.new_ab_add_moment(:,2:end)];
        Knee_combine.Total.keeogo.flexion_extension.values = [Knee_combine.Total.keeogo.flexion_extension.values, flex_ext_keeogo_tot];
        Knee_combine.Total.keeogo.ab_adduction.values = [Knee_combine.Total.keeogo.ab_adduction.values, ab_ad_duction_keeogo_tot];
    end
end
%%
a = fieldnames(Knee_combine.Right.control);
b = fieldnames(Knee_combine.Right);
% Calculate mean and standard deviation values
for i = 1:length(b)
    for j = 1:length(a)
        value_r = Knee_combine.Right.(b{i}).(a{j}).values;
        value_l = Knee_combine.Left.(b{i}).(a{j}).values;
        value_tot = Knee_combine.Total.(b{i}).(a{j}).values;
        Knee_combine.Right.(b{i}).(a{j}).mean = mean(value_r(:,1:end),2,"omitnan");
        Knee_combine.Left.(b{i}).(a{j}).mean = mean(value_l(:,1:end),2,"omitnan");
        Knee_combine.Total.(b{i}).(a{j}).mean = mean(value_tot(:,1:end),2,"omitnan");
        Knee_combine.Right.(b{i}).(a{j}).std = std(value_r(:,1:end),0,2,"omitnan");
        Knee_combine.Left.(b{i}).(a{j}).std = std(value_l(:,1:end),0,2,"omitnan");
        Knee_combine.Total.(b{i}).(a{j}).std = std(value_tot(:,1:end),0,2,"omitnan");
    end
end
%% plot
for i = 1:length(b)
    for j = 1:length(a)
        % Create +/- 1 StDev plot
        Knee_combine.Right.(b{i}).(a{j}).flex_curve1 = Knee_combine.Right.(b{i}).(a{j}).mean + Knee_combine.Right.(b{i}).(a{j}).std;
        Knee_combine.Left.(b{i}).(a{j}).flex_curve1 = Knee_combine.Left.(b{i}).(a{j}).mean + Knee_combine.Left.(b{i}).(a{j}).std;
        Knee_combine.Total.(b{i}).(a{j}).flex_curve1 = Knee_combine.Total.(b{i}).(a{j}).mean + Knee_combine.Total.(b{i}).(a{j}).std;
        Knee_combine.Right.(b{i}).(a{j}).flex_curve2 = Knee_combine.Right.(b{i}).(a{j}).mean - Knee_combine.Right.(b{i}).(a{j}).std;
        Knee_combine.Left.(b{i}).(a{j}).flex_curve2 = Knee_combine.Left.(b{i}).(a{j}).mean - Knee_combine.Left.(b{i}).(a{j}).std;
        Knee_combine.Total.(b{i}).(a{j}).flex_curve2 = Knee_combine.Total.(b{i}).(a{j}).mean - Knee_combine.Total.(b{i}).(a{j}).std;
        % Plot across x axis
        Knee_combine.Right.(b{i}).(a{j}).x = 1:numel(Knee_combine.Right.(b{i}).(a{j}).mean);
        Knee_combine.Left.(b{i}).(a{j}).x = 1:numel(Knee_combine.Left.(b{i}).(a{j}).mean);
        Knee_combine.Total.(b{i}).(a{j}).x = 1:numel(Knee_combine.Total.(b{i}).(a{j}).mean);
        % Fill in space between mean and StDev lines
        Knee_combine.Right.(b{i}).(a{j}).x2 = [Knee_combine.Right.(b{i}).(a{j}).x, fliplr(Knee_combine.Right.(b{i}).(a{j}).x)];
        Knee_combine.Left.(b{i}).(a{j}).x2 = [Knee_combine.Left.(b{i}).(a{j}).x, fliplr(Knee_combine.Left.(b{i}).(a{j}).x)];
        Knee_combine.Total.(b{i}).(a{j}).x2 = [Knee_combine.Total.(b{i}).(a{j}).x, fliplr(Knee_combine.Total.(b{i}).(a{j}).x)];
        Knee_combine.Right.(b{i}).(a{j}).inbetween = [Knee_combine.Right.(b{i}).(a{j}).flex_curve1', fliplr(Knee_combine.Right.(b{i}).(a{j}).flex_curve2')];
        Knee_combine.Left.(b{i}).(a{j}).inbetween = [Knee_combine.Left.(b{i}).(a{j}).flex_curve1', fliplr(Knee_combine.Left.(b{i}).(a{j}).flex_curve2')];
        Knee_combine.Total.(b{i}).(a{j}).inbetween = [Knee_combine.Total.(b{i}).(a{j}).flex_curve1', fliplr(Knee_combine.Total.(b{i}).(a{j}).flex_curve2')];
    end
end
% Flexion/Extension left
figure(1001)
hold on
xlim([0,101])
fill(Knee_combine.Left.control.flexion_extension.x2, -Knee_combine.Left.control.flexion_extension.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Left.control.flexion_extension.x, -Knee_combine.Left.control.flexion_extension.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Left.keeogo.flexion_extension.x2, -Knee_combine.Left.keeogo.flexion_extension.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Left.keeogo.flexion_extension.x, -Knee_combine.Left.keeogo.flexion_extension.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Left Flexion/Extension Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off

% Flexion/Extension right
figure(1002)
hold on
xlim([0,101])
fill(Knee_combine.Right.control.flexion_extension.x2, -Knee_combine.Right.control.flexion_extension.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Right.control.flexion_extension.x,-Knee_combine.Right.control.flexion_extension.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Right.keeogo.flexion_extension.x2, -Knee_combine.Right.keeogo.flexion_extension.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Right.keeogo.flexion_extension.x, -Knee_combine.Right.keeogo.flexion_extension.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Right Flexion/Extension Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off

% Abduction/Adduction left
figure(1003)
hold on
xlim([0,101])
fill(Knee_combine.Left.control.ab_adduction.x2, Knee_combine.Left.control.ab_adduction.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Left.control.ab_adduction.x, Knee_combine.Left.control.ab_adduction.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Left.keeogo.ab_adduction.x2, Knee_combine.Left.keeogo.ab_adduction.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Left.keeogo.ab_adduction.x, Knee_combine.Left.keeogo.ab_adduction.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Left Ab/Ad-duction Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off

% Abduction/Adduction right
figure(1004)
hold on
xlim([0,101])
fill(Knee_combine.Right.control.ab_adduction.x2, -Knee_combine.Right.control.ab_adduction.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Right.control.ab_adduction.x, -Knee_combine.Right.control.ab_adduction.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Right.keeogo.ab_adduction.x2, -Knee_combine.Right.keeogo.ab_adduction.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Right.keeogo.ab_adduction.x, -Knee_combine.Right.keeogo.ab_adduction.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Right Ab/Ad-duction Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off

% Flexion/Extension total
figure(1005)
hold on
xlim([0,101])
fill(Knee_combine.Total.control.flexion_extension.x2, -Knee_combine.Total.control.flexion_extension.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Total.control.flexion_extension.x,-Knee_combine.Total.control.flexion_extension.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Total.keeogo.flexion_extension.x2, -Knee_combine.Total.keeogo.flexion_extension.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Total.keeogo.flexion_extension.x, -Knee_combine.Total.keeogo.flexion_extension.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Total Flexion/Extension Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off

% Abduction/Adduction total
figure(1006)
hold on
xlim([0,101])
fill(Knee_combine.Total.control.ab_adduction.x2, Knee_combine.Total.control.ab_adduction.inbetween, 'b','FaceAlpha',0.4);
plot(Knee_combine.Total.control.ab_adduction.x, Knee_combine.Total.control.ab_adduction.mean, 'b', 'LineWidth', 2);
hold on
fill(Knee_combine.Total.keeogo.ab_adduction.x2, Knee_combine.Total.keeogo.ab_adduction.inbetween, 'r','FaceAlpha',0.4);
plot(Knee_combine.Total.keeogo.ab_adduction.x, Knee_combine.Total.keeogo.ab_adduction.mean, 'r', 'LineWidth', 2);
legend('Control StDev','Control Mean','Keeogo StDev','Keeogo Mean')
title('Total Ab/Ad-duction Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
hold off