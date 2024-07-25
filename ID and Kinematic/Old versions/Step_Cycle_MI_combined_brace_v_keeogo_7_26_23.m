% clc
% clearvars
function [RKnee, LKnee] = Step_Cycle_MI_combined_brace_v_keeogo_7_26_23(forceStruct,LKnee,RKnee,height,mass,trial_name,thresholds)
time = forceStruct.time;
ForceLT = abs(forceStruct.f1(:,2));
ForceRT = abs(forceStruct.f2(:,2));

zcr_threshold = 50; % zero crossin detection threshold / your value here

%% ForceLT processing
[~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds); % positive (pos) and negative (neg) slope crossing points
% nn_l = min(numel(heel_strike_data_l(:,1)),numel(toe_off_data_l(:,1)));
heel_time_l = cell2mat(heel_strike_data_l(:,1));
toe_time_l = cell2mat(toe_off_data_l(:,1));
nn_l = min(numel(heel_time_l(:,1)),numel(toe_time_l(:,1)));
for time_index = 1:nn_l
    if heel_time_l(time_index) > toe_time_l(time_index)
        for ci = 2:nn_l
            % best solution is to interpolate  the data on new time axis; 
            new_t1_l(1:100,ci) = linspace(heel_time_l(ci-1),toe_time_l(ci),100);
            LKnee.new_ForceLT(:,ci) = interp1(time,ForceLT,new_t1_l(:,ci)');
            LKnee.new_global_flexionmoment_l(:,ci) = interp1(time,LKnee.Moment(:,1),new_t1_l(:,ci)');
            LKnee.new_flexionmoment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,1),new_t1_l(:,ci)');
            LKnee.new_rotationmoment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,2),new_t1_l(:,ci)');
            LKnee.new_ab_add_moment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,3),new_t1_l(:,ci)');
            LKnee.new_ml_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,1),new_t1_l(:,ci)');
            LKnee.new_vertical_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,2),new_t1_l(:,ci)');
            LKnee.new_ap_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,3),new_t1_l(:,ci)');
            % new_t1_l = new_t1_l - new_t1_l(1); % put time axis start at zero
        end
    else
        for ci = 1:nn_l
            % best solution is to interpolate  the data on new time axis; 
            new_t1_l(1:100,ci) = linspace(heel_time_l(ci),toe_time_l(ci),100);
            LKnee.new_ForceLT(:,ci) = interp1(time,ForceLT,new_t1_l(:,ci)');
            LKnee.new_global_flexionmoment_l(:,ci) = interp1(time,LKnee.Moment(:,1),new_t1_l(:,ci)');
            LKnee.new_flexionmoment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,1),new_t1_l(:,ci)');
            LKnee.new_rotationmoment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,2),new_t1_l(:,ci)');
            LKnee.new_ab_add_moment_l(:,ci) = interp1(time,LKnee.left_seg_kneemoment(:,3),new_t1_l(:,ci)');
            LKnee.new_ml_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,1),new_t1_l(:,ci)');
            LKnee.new_vertical_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,2),new_t1_l(:,ci)');
            LKnee.new_ap_force_l(:,ci) = interp1(time,LKnee.left_seg_kneeforce(:,3),new_t1_l(:,ci)');
            % new_t1_l = new_t1_l - new_t1_l(1); % put time axis start at zero
        end
    end
end
a = fieldnames(LKnee);
for i = 4
    LKnee.(a{i}) = 100*LKnee.(a{i})/(mass*9.81*height);
end
normtime = linspace(0,100,100);
new_ForceLT_mean = mean(LKnee.new_ForceLT(:,2:end),2);
%
new_global_flexionmoment_l_mean = mean(LKnee.new_global_flexionmoment_l(:,2:end),2,"omitnan");
new_global_flexionmoment_l_std = std(LKnee.new_global_flexionmoment_l(:,2:end),0,2,"omitnan");
% figure(1);
% plot(normtime,new_ForceLT(:,2:end),normtime,new_ForceLT_mean,'*-k');
% title('ForceLT');
new_flex_moment_mean_l = mean(LKnee.new_flexionmoment_l(:,2:end),2,"omitnan");
new_flex_moment_std_l = std(LKnee.new_flexionmoment_l(:,2:end),0,2,"omitnan");
% figure(2);
% plot(normtime,-new_flexionmoment_l(:,2:end)/(height*mass*9.8)*100,normtime,-new_flex_moment_mean_l/(height*mass*9.8)*100,'*-k');
% title('Left Flexion Moment');xlabel('% Stance');ylabel('ML Moment (% HTxBW)')
new_rot_moment_mean_l = mean(LKnee.new_rotationmoment_l(:,2:end),2,"omitnan");
new_rot_moment_std_l = std(LKnee.new_rotationmoment_l(:,2:end),0,2,"omitnan");
% figure(3);
% plot(normtime,-new_rotationmoment_l(:,2:end)/(height*mass*9.8)*100,normtime,-new_rot_moment_mean_l/(height*mass*9.8)*100,'*-k');
% title('Left Rotation Moment');xlabel('% Stance');ylabel('Vertical Moment (% HTxBW)')
new_ab_add_moment_mean_l = mean(LKnee.new_ab_add_moment_l(:,2:end),2,"omitnan");
new_ab_add_moment_std_l = std(LKnee.new_ab_add_moment_l(:,2:end),0,2,"omitnan");
% figure(4);
% plot(normtime,new_ab_add_moment_l(:,2:end)/(height*mass*9.8)*100,normtime,new_ab_add_moment_mean_l/(height*mass*9.8)*100,'*-k');
% title('Left Ab/Adduction Moment');xlabel('% Stance');ylabel('AP Moment (% HTxBW)')
% new_ml_mean_l = mean(new_ml_force_l(:,2:end),2);
% figure(5);
% plot(normtime,new_ml_force_l(:,2:end)/(mass*9.8)*100,normtime,new_ml_mean_l/(mass*9.8)*100,'*-k');
% title('Left ML Force');xlabel('% Stance');ylabel('Force (% BW)')
% new_vertical_mean_l = mean(new_vertical_force_l(:,2:end),2);
% figure(6);
% plot(normtime,new_vertical_force_l(:,2:end)/(mass*9.8)*100,normtime,new_vertical_mean_l/(mass*9.8)*100,'*-k');
% title('Left Vertical Force');xlabel('% Stance');ylabel('Force (% BW)')
% new_ap_mean_l = mean(new_ap_force_l(:,2:end),2);
% figure(7);
% plot(normtime,-new_ap_force_l(:,2:end)/(mass*9.8)*100,normtime,-new_ap_mean_l/(mass*9.8)*100,'*-k');
% title('Left AP Force');xlabel('% Stance');ylabel('Force (% BW)')

%% ForceRT processing
[~, heel_strike_data,toe_off_data] = detect_cycles(ForceRT,time,thresholds);
heel_time_r = cell2mat(heel_strike_data(:,1));
toe_time_r = cell2mat(toe_off_data(:,1));
nn_r = min(numel(heel_time_r),numel(toe_time_r));
for time_index = 1:nn_r
    if heel_time_r(time_index) > toe_time_r(time_index)
        for ci = 2:nn_r
    % best solution is to interpolate  the data on new time axis; 
            new_t1_r(1:100,ci) = linspace(heel_time_r(ci-1),toe_time_r(ci),100)';
            RKnee.new_ForceRT(:,ci) = interp1(time,ForceRT,new_t1_r(:,ci)');
            RKnee.new_global_flexionmoment_r(:,ci) = interp1(time,RKnee.Moment(:,1),new_t1_r(:,ci)');
            RKnee.new_flexionmoment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,1),new_t1_r(:,ci)');
            RKnee.new_rotationmoment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,2),new_t1_r(:,ci)');
            RKnee.new_ab_add_moment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,3),new_t1_r(:,ci)');
            RKnee.new_ml_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,1),new_t1_r(:,ci)');
            RKnee.new_vertical_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,2),new_t1_r(:,ci)');
            RKnee.new_ap_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,3),new_t1_r(:,ci)');
            % new_t12 = (new_t1 - new_t1(1,ci))/new_t1(end,ci); % put time axis start at zero
            
        end
    else
        for ci = 1:nn_r
            new_t1_r(1:100,ci) = linspace(heel_time_r(ci),toe_time_r(ci),100)';
            RKnee.new_ForceRT(:,ci) = interp1(time,ForceRT,new_t1_r(:,ci)');
            RKnee.new_global_flexionmoment_r(:,ci) = interp1(time,RKnee.Moment(:,1),new_t1_r(:,ci)');
            RKnee.new_flexionmoment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,1),new_t1_r(:,ci)');
            RKnee.new_rotationmoment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,2),new_t1_r(:,ci)');
            RKnee.new_ab_add_moment_r(:,ci) = interp1(time,RKnee.right_seg_kneemoment(:,3),new_t1_r(:,ci)');
            RKnee.new_ml_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,1),new_t1_r(:,ci)');
            RKnee.new_vertical_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,2),new_t1_r(:,ci)');
            RKnee.new_ap_force_r(:,ci) = interp1(time,RKnee.right_seg_kneeforce(:,3),new_t1_r(:,ci)');
            % new_t12 = (n
        end
    end
end
a = fieldnames(RKnee);
for i = 4
    RKnee.(a{i}) = 100*RKnee.(a{i})/(mass*9.81*height);
end
% new_ForceRT_mean = mean(RKnee.new_ForceRT(:,2:end),2);
% %
% new_global_flexionmoment_mean = nanmean(RKnee.new_global_flexionmoment(:,2:end),2);
% new_global_flexionmoment_std = nanstd(RKnee.new_global_flexionmoment(:,2:end),0,2);
% % figure(14);
% % plot(normtime,new_ForceRT(:,2:end),normtime,new_ForceRT_mean,'*-k');
% % title('ForceRT');
% new_flex_moment_mean = nanmean(RKnee.new_flexionmoment(:,2:end),2);
% new_flex_moment_std = nanstd(RKnee.new_flexionmoment(:,2:end),0,2);
% % figure(15);
% % plot(normtime,-new_flexionmoment(:,2:end)/(height*mass*9.8)*100,normtime,-new_flex_moment_mean/(height*mass*9.8)*100,'*-k');
% % title('Right Flexion Moment');xlabel('% Stance');ylabel('ML Moment (% HTxBW)')
% new_rot_moment_mean = nanmean(RKnee.new_rotationmoment(:,2:end),2);
% new_rot_moment_std = nanstd(RKnee.new_rotationmoment(:,2:end),0,2);
% % figure(16);
% % plot(normtime,-new_rotationmoment(:,2:end)/(height*mass*9.8)*100,normtime,-new_rot_moment_mean/(height*mass*9.8)*100,'*-k');
% % title('Right Rotation Moment');xlabel('% Stance');ylabel('Vertical Moment (% HTxBW)')
% new_ab_add_moment_mean = nanmean(RKnee.new_ab_add_moment(:,2:end),2);
% new_ab_add_moment_std = nanstd(RKnee.new_ab_add_moment(:,2:end),0,2);
% % figure(17);
% % plot(normtime,new_ab_add_moment(:,2:end)/(height*mass*9.8)*100,normtime,new_ab_add_moment_mean/(height*mass*9.8)*100,'*-k');
% % title('Right Ab/Adduction Moment');xlabel('% Stance');ylabel('AP Moment (% HTxBW)')
% % new_ml_mean = mean(new_ml_force(:,2:end),2);
% % figure(18);
% % plot(normtime,new_ml_force(:,2:end)/(mass*9.8)*100,normtime,new_ml_mean/(mass*9.8)*100,'*-k');
% % title('Right ML Force');xlabel('% Stance');ylabel('Force (% BW)')
% % new_vertical_mean = mean(new_vertical_force(:,2:end),2);
% % figure(19);
% % plot(normtime,new_vertical_force(:,2:end)/(mass*9.8)*100,normtime,new_vertical_mean/(mass*9.8)*100,'*-k');
% % title('Right Vertical Force');xlabel('% Stance');ylabel('Force (% BW)')
% % new_ap_mean = mean(new_ap_force(:,2:end),2);
% % figure(20);
% % plot(normtime,-new_ap_force(:,2:end)/(mass*9.8)*100,normtime,-new_ap_mean/(mass*9.8)*100,'*-k');
% % title('Right AP Force');xlabel('% Stance');ylabel('Force (% BW)')
% 
% %% Standard deviation plot
% %make shorter
% figure(100)
% hold on
% flex_curve1 = -new_flex_moment_mean*100/(height*mass*9.8) + new_flex_moment_std*100/(height*mass*9.8);
% flex_curve2 = -new_flex_moment_mean*100/(height*mass*9.8) - new_flex_moment_std*100/(height*mass*9.8);
% x = 1:numel(new_flex_moment_mean);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, -new_flex_moment_mean*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, -new_flex_moment_mean*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on
% end
% title('Flexion/Extension Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)');
% hold off
% 
% figure(101)
% hold on
% flex_curve1 = -new_flex_moment_mean_l*100/(height*mass*9.8) + new_flex_moment_std_l*100/(height*mass*9.8);
% flex_curve2 = -new_flex_moment_mean_l*100/(height*mass*9.8) - new_flex_moment_std_l*100/(height*mass*9.8);
% x = 1:numel(new_flex_moment_mean_l);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% hold on;
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, -new_flex_moment_mean_l*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, -new_flex_moment_mean_l*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
% end
% title('Left Flexion/Extension Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(102)
% hold on
% flex_curve1 = -new_rot_moment_mean*100/(height*mass*9.8) + new_rot_moment_std*100/(height*mass*9.8);
% flex_curve2 = -new_rot_moment_mean*100/(height*mass*9.8) - new_rot_moment_std*100/(height*mass*9.8);
% x = 1:numel(new_rot_moment_mean);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, -new_rot_moment_mean*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on;
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, -new_rot_moment_mean*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on;
% end
% title('Rotational Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(103)
% hold on
% flex_curve1 = -new_rot_moment_mean_l*100/(height*mass*9.8) + new_rot_moment_std_l*100/(height*mass*9.8);
% flex_curve2 = -new_rot_moment_mean_l*100/(height*mass*9.8) - new_rot_moment_std_l*100/(height*mass*9.8);
% x = 1:numel(new_rot_moment_mean_l);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, -new_rot_moment_mean_l*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on;
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, -new_rot_moment_mean_l*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on;
% end
% title('Left Rotational Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(104)
% hold on
% flex_curve1 = new_ab_add_moment_mean*100/(height*mass*9.8) + new_ab_add_moment_std*100/(height*mass*9.8);
% flex_curve2 = new_ab_add_moment_mean*100/(height*mass*9.8) - new_ab_add_moment_std*100/(height*mass*9.8);
% x = 1:numel(new_ab_add_moment_mean);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, new_ab_add_moment_mean*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, new_ab_add_moment_mean*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on
% end
% title('Ab/Adduction Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(105)
% hold on
% flex_curve1 = new_ab_add_moment_mean_l*100/(height*mass*9.8) + new_ab_add_moment_std_l*100/(height*mass*9.8);
% flex_curve2 = new_ab_add_moment_mean_l*100/(height*mass*9.8) - new_ab_add_moment_std_l*100/(height*mass*9.8);
% x = 1:numel(new_ab_add_moment_mean_l);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, new_ab_add_moment_mean_l*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Neutral StDev','Neutral mean')
%     hold on
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, new_ab_add_moment_mean_l*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Uphill StDev','Uphill Mean')
%     hold on
% end
% title('Left Ab/Adduction Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(106)
% hold on
% flex_curve1 = new_global_flexionmoment_mean*100/(height*mass*9.8) + new_global_flexionmoment_std*100/(height*mass*9.8);
% flex_curve2 = new_global_flexionmoment_mean*100/(height*mass*9.8) - new_global_flexionmoment_std*100/(height*mass*9.8);
% x = 1:numel(new_global_flexionmoment_mean);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, new_global_flexionmoment_mean*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, new_global_flexionmoment_mean*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on
% 
% end
% title('Right Global Flexion Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
% 
% figure(107)
% hold on
% flex_curve1 = new_global_flexionmoment_l_mean*100/(height*mass*9.8) + new_global_flexionmoment_l_std*100/(height*mass*9.8);
% flex_curve2 = new_global_flexionmoment_l_mean*100/(height*mass*9.8) - new_global_flexionmoment_l_std*100/(height*mass*9.8);
% x = 1:numel(new_global_flexionmoment_mean);
% x2 = [x, fliplr(x)];
% inBetween = [flex_curve1', fliplr(flex_curve2')];
% 
% if contains(lower(trial_name),'keeogo')
%     fill(x2, inBetween, 'r','FaceAlpha',0.4);
%     plot(x, new_global_flexionmoment_l_mean*100/(height*mass*9.8), 'r', 'LineWidth', 2);
%     legend('Keeogo StDev','Keeogo mean')
%     hold on
% else
%     fill(x2, inBetween, 'b','FaceAlpha',0.4);
%     plot(x, new_global_flexionmoment_l_mean*100/(height*mass*9.8), 'b', 'LineWidth', 2);
%     legend('Control StDev','Control Mean')
%     hold on
% end
% title('Left Global Flexion Moment');xlabel('% Stance');ylabel('Moment (%BW*HT)')
% hold off
end
