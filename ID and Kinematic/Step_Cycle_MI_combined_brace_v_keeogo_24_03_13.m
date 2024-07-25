% clc
% clearvars
function [Knee] = Step_Cycle_MI_combined_brace_v_keeogo_24_03_13(forceStruct,Knee,height,mass,thresholds)
time = forceStruct.time;
ForceLT = abs(forceStruct.f1(:,2));
ForceRT = abs(forceStruct.f2(:,2));

%% Left Knee processing
% Determine individual cycles
[~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds); % positive (pos) and negative (neg) slope crossing points

% Convert data from cell format to matrix
heel_time_l = cell2mat(heel_strike_data_l(:,1));
toe_time_l = cell2mat(toe_off_data_l(:,1));
% min # number of steps detected
nn_l = min(numel(heel_time_l(:,1)),numel(toe_time_l(:,1))); 

new_t1_l = zeros(101,nn_l);
% Record normalized parameter values per stance
for i = 1:nn_l
    if heel_time_l(i) > toe_time_l(i)
        for ci = 2:nn_l
            new_t1_l(1:101,ci) = linspace(heel_time_l(ci-1),toe_time_l(ci),101);
        end
    else
        for ci = 1:nn_l
            new_t1_l(1:101,ci) = linspace(heel_time_l(ci),toe_time_l(ci),101);
        end
    end
end

for k = 1:nn_l
    % best solution is to interpolate  the data on new time axis; 
    Knee.Left.new_ForceLT(:,k) = interp1(time,ForceLT,new_t1_l(:,k)');
    Knee.Left.new_global_flexionmoment(:,k) = interp1(time,Knee.Left.Moment(:,1),new_t1_l(:,k)');
    Knee.Left.new_flexionmoment(:,k) = interp1(time,Knee.Left.seg_kneemoment(:,1),new_t1_l(:,k)');
    Knee.Left.new_rotationmoment(:,k) = interp1(time,Knee.Left.seg_kneemoment(:,2),new_t1_l(:,k)');
    Knee.Left.new_ab_add_moment(:,k) = interp1(time,Knee.Left.seg_kneemoment(:,3),new_t1_l(:,k)');
    Knee.Left.new_ml_force(:,k) = interp1(time,Knee.Left.seg_kneeforce(:,1),new_t1_l(:,k)');
    Knee.Left.new_vertical_force(:,k) = interp1(time,Knee.Left.seg_kneeforce(:,2),new_t1_l(:,k)');
    Knee.Left.new_ap_force(:,k) = interp1(time,Knee.Left.seg_kneeforce(:,3),new_t1_l(:,k)');
end


parameters = {'new_global_flexionmoment','new_flexionmoment',...
    'new_rotationmoment','new_ab_add_moment'};
force_param = {'new_ForceLT','new_ml_force','new_vertical_force',...
    'new_ap_force'};
% % Normalize to BW*Ht
for i = 1:length(parameters)
    Knee.Left.(parameters{i}) = 100*Knee.Left.(parameters{i})/(mass*9.81*height);
end
for i = 1:length(parameters)
    Knee.Left.(force_param{i}) = 100*Knee.Left.(force_param{i})/(mass*9.81);
end

%% Right Knee processing
[~, heel_strike_data,toe_off_data] = detect_cycles(ForceRT,time,thresholds);
heel_time_r = cell2mat(heel_strike_data(:,1));
toe_time_r = cell2mat(toe_off_data(:,1));

% min # number of steps detected
nn_r = min(numel(heel_time_r),numel(toe_time_r));
new_t1_r = zeros(101, nn_r);
% Record normalized parameter values per stance
for i = 1:nn_r
    if heel_time_r(i) > toe_time_r(i)
        for ci = 2:nn_r
            new_t1_r(1:101,ci) = linspace(heel_time_r(ci-1),toe_time_r(ci),101)';
        end
    else
        for ci = 1:nn_r
            new_t1_r(1:101,ci) = linspace(heel_time_r(ci),toe_time_r(ci),101)';
        end
    end
end

for k = 1:nn_r
    Knee.Right.new_ForceRT(:,k) = interp1(time,ForceRT,new_t1_r(:,k)');
    Knee.Right.new_global_flexionmoment(:,k) = interp1(time,Knee.Right.Moment(:,1),new_t1_r(:,k)');
    Knee.Right.new_flexionmoment(:,k) = interp1(time,Knee.Right.seg_kneemoment(:,1),new_t1_r(:,k)');
    Knee.Right.new_rotationmoment(:,k) = interp1(time,Knee.Right.seg_kneemoment(:,2),new_t1_r(:,k)');
    Knee.Right.new_ab_add_moment(:,k) = interp1(time,Knee.Right.seg_kneemoment(:,3),new_t1_r(:,k)');
    Knee.Right.new_ml_force(:,k) = interp1(time,Knee.Right.seg_kneeforce(:,1),new_t1_r(:,k)');
    Knee.Right.new_vertical_force(:,k) = interp1(time,Knee.Right.seg_kneeforce(:,2),new_t1_r(:,k)');
    Knee.Right.new_ap_force(:,k) = interp1(time,Knee.Right.seg_kneeforce(:,3),new_t1_r(:,k)');
end

parameters = {'new_global_flexionmoment','new_flexionmoment',...
    'new_rotationmoment','new_ab_add_moment'};
force_param = {'new_ForceRT','new_ml_force','new_vertical_force',...
    'new_ap_force'};
% Normalize by BW*Ht
for i = 1:length(parameters)
    Knee.Right.(parameters{i}) = 100*Knee.Right.(parameters{i})/(mass*9.81*height);
end
% Normalize by BW
for i = 1:length(force_param)
    Knee.Right.(force_param{i}) = 100*Knee.Right.(force_param{i})/(mass*9.81);
end

end
