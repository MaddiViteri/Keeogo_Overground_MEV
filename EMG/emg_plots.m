%% plot EMG data
%%
function [stuff] = emg_plots(gait_emg_metrics_combined1_rms,gait_emg_metrics_combined1,gait_emg_swing_metrics_combined)
%%
X = categorical({'Loading','Early','Mid','Late'});
X = reordercats(X,{'Loading','Early','Mid','Late'});
Z = categorical({'Quad','Hamstring','Gas'});
Z = reordercats(Z,{'Quad','Hamstring','Gas'});
%% CCI
CCI_mean_right = [gait_emg_metrics_combined1_rms.control.CCI_r_loading_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_r_loading_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_r_early_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_r_early_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_r_mid_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_r_mid_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_r_late_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_r_late_rms.mean];

CCI_mean_left = [gait_emg_metrics_combined1_rms.control.CCI_l_loading_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_l_loading_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_l_early_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_l_early_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_l_mid_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_l_mid_rms.mean;...
    gait_emg_metrics_combined1_rms.control.CCI_l_late_rms.mean, gait_emg_metrics_combined1_rms.keeogo.CCI_l_late_rms.mean];

CCI_std_right = [gait_emg_metrics_combined1_rms.control.CCI_r_loading_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_r_loading_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_r_early_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_r_early_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_r_mid_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_r_mid_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_r_late_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_r_late_rms.std];

CCI_std_left = [gait_emg_metrics_combined1_rms.control.CCI_l_loading_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_l_loading_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_l_early_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_l_early_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_l_mid_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_l_mid_rms.std;...
    gait_emg_metrics_combined1_rms.control.CCI_l_late_rms.std, gait_emg_metrics_combined1_rms.keeogo.CCI_l_late_rms.std];
%% CCR
CCR_zero_right = [gait_emg_metrics_combined1.control.CCR_zero_r_loading.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_r_loading.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_r_early.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_r_early.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_r_mid.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_r_mid.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_r_late.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_r_late.mean];

CCR_zero_left = [gait_emg_metrics_combined1.control.CCR_zero_l_loading.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_l_loading.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_l_early.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_l_early.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_l_mid.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_l_mid.mean;...
    gait_emg_metrics_combined1.control.CCR_zero_l_late.mean, gait_emg_metrics_combined1.keeogo.CCR_zero_l_late.mean];

CCR_zero_std_right = [gait_emg_metrics_combined1.control.CCR_zero_r_loading.std, gait_emg_metrics_combined1.keeogo.CCR_zero_r_loading.std;...
    gait_emg_metrics_combined1.control.CCR_zero_r_early.std, gait_emg_metrics_combined1.keeogo.CCR_zero_r_early.std;...
    gait_emg_metrics_combined1.control.CCR_zero_r_mid.std, gait_emg_metrics_combined1.keeogo.CCR_zero_r_mid.std;...
    gait_emg_metrics_combined1.control.CCR_zero_r_late.std, gait_emg_metrics_combined1.keeogo.CCR_zero_r_late.std];

CCR_zero_std_left = [gait_emg_metrics_combined1.control.CCR_zero_l_loading.std, gait_emg_metrics_combined1.keeogo.CCR_zero_l_loading.std;...
    gait_emg_metrics_combined1.control.CCR_zero_l_early.std, gait_emg_metrics_combined1.keeogo.CCR_zero_l_early.std;...
    gait_emg_metrics_combined1.control.CCR_zero_l_mid.std, gait_emg_metrics_combined1.keeogo.CCR_zero_l_mid.std;...
    gait_emg_metrics_combined1.control.CCR_zero_l_late.std, gait_emg_metrics_combined1.keeogo.CCR_zero_l_late.std];
%% Swing Right
swing_full = [gait_emg_swing_metrics_combined.control.rquad_swing.mean, gait_emg_swing_metrics_combined.keeogo.rquad_swing.mean;...
    gait_emg_swing_metrics_combined.control.rham_swing.mean, gait_emg_swing_metrics_combined.keeogo.rham_swing.mean;...
    gait_emg_swing_metrics_combined.control.rgas_swing.mean, gait_emg_swing_metrics_combined.keeogo.rgas_swing.mean];
swing_full_std = [gait_emg_swing_metrics_combined.control.rquad_swing.std, gait_emg_swing_metrics_combined.keeogo.rquad_swing.std;...
    gait_emg_swing_metrics_combined.control.rham_swing.std, gait_emg_swing_metrics_combined.keeogo.rham_swing.std;...
    gait_emg_swing_metrics_combined.control.rgas_swing.std, gait_emg_swing_metrics_combined.keeogo.rgas_swing.std];
% Prep
swing_prep = [gait_emg_swing_metrics_combined.control.rquad_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.rquad_stride_prep.mean;...
    gait_emg_swing_metrics_combined.control.rham_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.rham_stride_prep.mean;...
    gait_emg_swing_metrics_combined.control.rgas_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.rgas_stride_prep.mean];
swing_prep_std = [gait_emg_swing_metrics_combined.control.rquad_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.rquad_stride_prep.std;...
    gait_emg_swing_metrics_combined.control.rham_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.rham_stride_prep.std;...
    gait_emg_swing_metrics_combined.control.rgas_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.rgas_stride_prep.std];
%% Swing left
swing_full_left = [gait_emg_swing_metrics_combined.control.lquad_swing.mean, gait_emg_swing_metrics_combined.keeogo.lquad_swing.mean;...
    gait_emg_swing_metrics_combined.control.lham_swing.mean, gait_emg_swing_metrics_combined.keeogo.lham_swing.mean;...
    gait_emg_swing_metrics_combined.control.lgas_swing.mean, gait_emg_swing_metrics_combined.keeogo.lgas_swing.mean];
swing_full_left_std = [gait_emg_swing_metrics_combined.control.lquad_swing.std, gait_emg_swing_metrics_combined.keeogo.lquad_swing.std;...
    gait_emg_swing_metrics_combined.control.lham_swing.std, gait_emg_swing_metrics_combined.keeogo.lham_swing.std;...
    gait_emg_swing_metrics_combined.control.lgas_swing.std, gait_emg_swing_metrics_combined.keeogo.lgas_swing.std];
%%Prep
swing_prep_left = [gait_emg_swing_metrics_combined.control.lquad_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.lquad_stride_prep.mean;...
    gait_emg_swing_metrics_combined.control.lham_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.lham_stride_prep.mean;...
    gait_emg_swing_metrics_combined.control.lgas_stride_prep.mean, gait_emg_swing_metrics_combined.keeogo.lgas_stride_prep.mean];
swing_prep_left_std = [gait_emg_swing_metrics_combined.control.lquad_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.lquad_stride_prep.std;...
    gait_emg_swing_metrics_combined.control.lham_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.lham_stride_prep.std;...
    gait_emg_swing_metrics_combined.control.lgas_stride_prep.std, gait_emg_swing_metrics_combined.keeogo.lgas_stride_prep.std];
%% Stance Sum
sum_left = [gait_emg_metrics_combined1_rms.control.sum_l_loading_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_l_loading_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_l_early_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_l_early_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_l_mid_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_l_mid_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_l_late_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_l_late_rms.mean];
sum_left_std = [gait_emg_metrics_combined1_rms.control.sum_l_loading_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_l_loading_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_l_early_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_l_early_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_l_mid_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_l_mid_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_l_late_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_l_late_rms.std];
sum_right = [gait_emg_metrics_combined1_rms.control.sum_r_loading_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_r_loading_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_r_early_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_r_early_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_r_mid_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_r_mid_rms.mean;...
    gait_emg_metrics_combined1_rms.control.sum_r_late_rms.mean, gait_emg_metrics_combined1_rms.keeogo.sum_r_late_rms.mean];
sum_right_std = [gait_emg_metrics_combined1_rms.control.sum_r_loading_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_r_loading_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_r_early_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_r_early_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_r_mid_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_r_mid_rms.std;...
    gait_emg_metrics_combined1_rms.control.sum_r_late_rms.std, gait_emg_metrics_combined1_rms.keeogo.sum_r_late_rms.std];
%% Plot CCI
figure;
b  = bar(X,CCI_mean_right);
legend('Control','Keeogo')
ylabel('CCI')
title('CCI Right')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(CCI_mean_right);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',CCI_mean_right,CCI_std_right,'k','linestyle','none');
hold off

figure;
b  = bar(X,CCI_mean_left);
legend('Control','Keeogo')
ylabel('CCI')
title('CCI Left')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(CCI_mean_left);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',CCI_mean_left,CCI_std_left,'k','linestyle','none');
hold off
stuff=[];
%% Plot CCR
figure;
b  = bar(X,CCR_zero_right);
legend('Control','Keeogo')
ylabel('CCR')
title('CCR Right')
ylim([-1, 1]);
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(CCR_zero_right);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',CCR_zero_right,CCR_zero_std_right,'k','linestyle','none');
hold off

figure;
b  = bar(X,CCR_zero_left);
legend('Control','Keeogo')
ylabel('CCR')
title('CCR Left')
ylim([-1, 1]);
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(CCR_zero_left);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',CCR_zero_left,CCR_zero_std_left,'k','linestyle','none');
hold off
stuff=[];
%% Plot Swing Right
figure;
b  = bar(Z,swing_full);
legend('Control','Keeogo')
ylabel('RMS (%max)')
title('Swing Full Right')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(swing_full);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',swing_full,swing_full_std,'k','linestyle','none');
hold off
%%Prep
figure;
b  = bar(Z,swing_prep);
legend('Control','Keeogo')
ylabel('RMS (%max)')
title('Swing Prep Right')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(swing_prep);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',swing_prep,swing_prep_std,'k','linestyle','none');
hold off
stuff=[];
%% Plot Swing Left
figure;
b  = bar(Z,swing_full_left);
legend('Control','Keeogo')
ylabel('RMS (%max)')
title('Swing Left Full')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(swing_full_left);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',swing_full_left,swing_full_left_std,'k','linestyle','none');
hold off

%%Prep
figure;
b  = bar(Z,swing_prep_left);
legend('Control','Keeogo')
ylabel('RMS (%max)')
title('Swing Left Prep')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(swing_prep_left);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',swing_prep_left,swing_prep_left_std,'k','linestyle','none');
hold off
stuff=[];
%% Plot Sum
figure;
b  = bar(X,sum_left);
legend('Control','Keeogo')
ylabel('Sum')
title('Sum Left')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(sum_left);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',sum_left,sum_left_std,'k','linestyle','none');
hold off

figure;
b  = bar(X,sum_right);
legend('Control','Keeogo')
ylabel('Sum')
title('Sum Right')
    % b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(sum_right);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',sum_right,sum_right_std,'k','linestyle','none');
hold off
stuff=[];