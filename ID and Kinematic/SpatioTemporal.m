% markerStruct = trial.neutral_control.markers;
% forceStruct = trial.neutral_control.force;
% trial_file = platform.neutral_control10001;
% thresholds = 50;
%%
function [gait_parameters] = SpatioTemporal(forceStruct,markerStruct,trial_file,thresholds)

% Rename variables
time = forceStruct.time;
ForceLT = forceStruct.f1(:,2);
ForceRT = forceStruct.f2(:,2); 
gait_parameters = struct();

%% Left
% Determine contacts
[~, heel_strike_data_l,toe_off_data_l] = detect_cycles(ForceLT,time,thresholds);
heel_time_l = round(cell2mat(heel_strike_data_l(:,1))*1000);
toe_time_l = round(cell2mat(toe_off_data_l(:,1))*1000);
% Upsample platform data to match force data length

trial_fileLength = height(trial_file);

disp(trial_fileLength);

ForceLTLength = height(ForceLT);

trial_file = table2array(trial_file);

class(trial_file); class(trial_fileLength); class(ForceLTLength);

test = ScaleTime(trial_file,1,trial_fileLength,ForceLTLength);
% min # number of steps detected
nn_l = min(numel(heel_time_l),numel(toe_time_l));

% Record normalized parameter values per stance
    k= 1;
for ci = 2:nn_l

    if heel_time_l(ci) > toe_time_l(ci)
        gait_parameters.left.stride_length(k) = markerStruct.LHEE(heel_time_l(ci),3)-markerStruct.LHEE(heel_time_l(ci-1),3)+mean([test(heel_time_l(ci),2),test(heel_time_l(ci-1),2)])*(test(heel_time_l(ci),1)-test(heel_time_l(ci-1),1));
        gait_parameters.left.stride_time(k) = test(heel_time_l(ci),1)-test(heel_time_l(ci-1),1);
        % best solution is to interpolate  the data on new time axis; 
        gait_parameters.left.stancetime(k) = time(toe_time_l(ci))-time(heel_time_l(ci-1));
        gait_parameters.left.swingtime(k) = time(heel_time_l(ci))-time(toe_time_l(ci));
        gait_parameters.left.stance_perc(k) = gait_parameters.left.stancetime(k)/(gait_parameters.left.stancetime(k)+gait_parameters.left.swingtime(k));
        gait_parameters.left.swing_perc(k) = gait_parameters.left.swingtime(k)/(gait_parameters.left.stancetime(k)+gait_parameters.left.swingtime(k));
        gait_parameters.left.toe_out(k) = 180/pi*mean(unwrap(atan2(markerStruct.L2Meta((heel_time_l(ci-1)+round(.15*(toe_time_l(ci)-heel_time_l(ci-1)))):(heel_time_l(ci-1)+round(.5*(toe_time_l(ci)-heel_time_l(ci-1)))),1)...
            -markerStruct.LHEE((heel_time_l(ci-1)+round(.15*(toe_time_l(ci)-heel_time_l(ci-1)))):(heel_time_l(ci-1)+round(.5*(toe_time_l(ci)-heel_time_l(ci-1)))),1),...
            markerStruct.L2Meta((heel_time_l(ci-1)+round(.15*(toe_time_l(ci)-heel_time_l(ci-1)))):(heel_time_l(ci-1)+round(.5*(toe_time_l(ci)-heel_time_l(ci-1)))),3)...
            -markerStruct.LHEE((heel_time_l(ci-1)+round(.15*(toe_time_l(ci)-heel_time_l(ci-1)))):(heel_time_l(ci-1)+round(.5*(toe_time_l(ci)-heel_time_l(ci-1)))),3))));
        k=k+1;
        
    else
        gait_parameters.left.stride_length(k) = markerStruct.LHEE(heel_time_l(ci),3)-markerStruct.LHEE(heel_time_l(ci-1),3)+mean([test(heel_time_l(ci),2),test(heel_time_l(ci-1),2)])*(test(heel_time_l(ci),1)-test(heel_time_l(ci-1),1));
        gait_parameters.left.stride_time(k) = test(heel_time_l(ci),1)-test(heel_time_l(ci-1),1);
        gait_parameters.left.stancetime(k) = time(toe_time_l(ci))-time(heel_time_l(ci));
        gait_parameters.left.swingtime(k) = time(heel_time_l(ci))-time(toe_time_l(ci-1));
        gait_parameters.left.stance_perc(k) = gait_parameters.left.stancetime(k)/(gait_parameters.left.stancetime(k)+gait_parameters.left.swingtime(k));
        gait_parameters.left.swing_perc(k) = gait_parameters.left.swingtime(k)/(gait_parameters.left.stancetime(k)+gait_parameters.left.swingtime(k));
        gait_parameters.left.toe_out(k) = 180/pi*mean(unwrap(atan2(markerStruct.L2Meta((heel_time_l(ci)+round(.15*(toe_time_l(ci)-heel_time_l(ci)))):(heel_time_l(ci)+round(.5*(toe_time_l(ci)-heel_time_l(ci)))),1)...
            -markerStruct.LHEE((heel_time_l(ci)+round(.15*(toe_time_l(ci)-heel_time_l(ci)))):(heel_time_l(ci)+round(.5*(toe_time_l(ci)-heel_time_l(ci)))),1),...
            markerStruct.L2Meta((heel_time_l(ci)+round(.15*(toe_time_l(ci)-heel_time_l(ci)))):(heel_time_l(ci)+round(.5*(toe_time_l(ci)-heel_time_l(ci)))),3)...
            -markerStruct.LHEE((heel_time_l(ci)+round(.15*(toe_time_l(ci)-heel_time_l(ci)))):(heel_time_l(ci)+round(.5*(toe_time_l(ci)-heel_time_l(ci)))),3))));
        k=k+1;

    end
end

%% Right
% Detect cycles
[~, heel_strike_data_r,toe_off_data_r] = detect_cycles(ForceRT,time,thresholds);
heel_time_r = round(cell2mat(heel_strike_data_r(:,1))*1000);
toe_time_r = round(cell2mat(toe_off_data_r(:,1))*1000);
nn_r = min(numel(heel_time_r),numel(toe_time_r));

% Record normalized parameter values per stance
k=1;
for ci = 2:nn_r
    
    if heel_time_r(ci) > toe_time_r(ci)
        
            gait_parameters.right.stride_length(k) = markerStruct.RHEE(heel_time_r(ci),3)-markerStruct.RHEE(heel_time_r(ci-1),3)+mean([test(heel_time_r(ci),4),test(heel_time_r(ci-1),4)])*(test(heel_time_r(ci),1)-test(heel_time_r(ci-1),1));
            gait_parameters.right.stride_time(k) = test(heel_time_r(ci),1)-test(heel_time_r(ci-1),1);
            % best solution is to interpolate  the data on new time axis; 
            gait_parameters.right.stancetime(k) = time(toe_time_r(ci))-time(heel_time_r(ci-1));
            gait_parameters.right.swingtime(k) = time(heel_time_r(ci))-time(toe_time_r(ci));
            gait_parameters.right.stance_perc(k) = gait_parameters.right.stancetime_l(k)/(gait_parameters.right.stancetime_l(k)+gait_parameters.right.swingtime_r(k));
            gait_parameters.right.swing_perc(k) = gait_parameters.swingtime_r(k)/(gait_parameters.right.stancetime_l(k)+gait_parameters.right.swingtime_r(k));
            gait_parameters.right.toe_out(k) = 180/pi*mean(unwrap(atan2(markerStruct.R2Meta((heel_time_r(ci-1)+round(.15*(toe_time_r(ci)-heel_time_r(ci-1)))):(heel_time_r(ci-1)+round(.5*(toe_time_r(ci)-heel_time_r(ci-1)))),1)...
                -markerStruct.RHEE((heel_time_r(ci-1)+round(.15*(toe_time_r(ci)-heel_time_r(ci-1)))):(heel_time_r(ci-1)+round(.5*(toe_time_r(ci)-heel_time_r(ci-1)))),1),...
                markerStruct.R2Meta((heel_time_r(ci-1)+round(.15*(toe_time_r(ci)-heel_time_r(ci-1)))):(heel_time_r(ci-1)+round(.5*(toe_time_r(ci)-heel_time_r(ci-1)))),3)...
                -markerStruct.RHEE((heel_time_r(ci-1)+round(.15*(toe_time_r(ci)-heel_time_r(ci-1)))):(heel_time_r(ci-1)+round(.5*(toe_time_r(ci)-heel_time_r(ci-1)))),3))));
            
            k=k+1;
    else

            gait_parameters.right.stride_length(k) = markerStruct.RHEE(heel_time_r(ci),3)...
                -markerStruct.RHEE(heel_time_r(ci-1),3)+...
                mean([test(heel_time_r(ci),4),test(heel_time_r(ci-1),4)])*(test(heel_time_r(ci),1)...
                -test(heel_time_r(ci-1),1));
            gait_parameters.right.stride_time(k) = test(heel_time_r(ci),1)-test(heel_time_r(ci-1),1);
            gait_parameters.right.stancetime(k) = time(toe_time_r(ci))-time(heel_time_r(ci));
            gait_parameters.right.swingtime(k) = time(heel_time_r(ci))-time(toe_time_r(ci-1));
            gait_parameters.right.stance_perc(k) = gait_parameters.right.stancetime(k)/(gait_parameters.right.stancetime(k)+gait_parameters.right.swingtime(k));
            gait_parameters.right.swing_perc(k) = gait_parameters.right.swingtime(k)/(gait_parameters.right.stancetime(k)+gait_parameters.right.swingtime(k));
            gait_parameters.right.toe_out(k) = 180/pi*mean(unwrap(atan2(markerStruct.R2Meta((heel_time_r(ci)+round(.15*(toe_time_r(ci)-heel_time_r(ci)))):(heel_time_r(ci)+round(.5*(toe_time_r(ci)-heel_time_r(ci)))),1)...
                -markerStruct.RHEE((heel_time_r(ci)+round(.15*(toe_time_r(ci)-heel_time_r(ci)))):(heel_time_r(ci)+round(.5*(toe_time_r(ci)-heel_time_r(ci)))),1),...
                markerStruct.R2Meta((heel_time_r(ci)+round(.15*(toe_time_r(ci)-heel_time_r(ci)))):(heel_time_r(ci)+round(.5*(toe_time_r(ci)-heel_time_r(ci)))),3)...
                -markerStruct.RHEE((heel_time_r(ci)+round(.15*(toe_time_r(ci)-heel_time_r(ci)))):(heel_time_r(ci)+round(.5*(toe_time_r(ci)-heel_time_r(ci)))),3))));
            k=k+1;
    end
end

%% Stride length combo
gait_parameters.total.stride_length= [gait_parameters.right.stride_length,gait_parameters.left.stride_length];
%% Stride time combine
gait_parameters.total.stride_time = [gait_parameters.right.stride_time, gait_parameters.left.stride_time];
%% Stance time combined
gait_parameters.total.stancetime = [gait_parameters.right.stancetime, gait_parameters.left.stancetime];
%% Swing time combined
gait_parameters.total.swingtime = [gait_parameters.right.swingtime, gait_parameters.left.swingtime];
%% Stance percentage combined
gait_parameters.total.stance_perc = [gait_parameters.right.stance_perc, gait_parameters.left.stance_perc];
%% Swing Percentage combined
gait_parameters.total.swing_perc = [gait_parameters.right.swing_perc, gait_parameters.left.swing_perc];
%% Toe out combined
gait_parameters.total.toe_out = [-(gait_parameters.right.toe_out),(gait_parameters.left.toe_out)];
%% Step Width and Length
k=1;
for i = 1:(length(heel_time_r)-1)
    gait_parameters.right.step_width(k) =  markerStruct.LHEE(heel_time_r(i),1)-markerStruct.RHEE(heel_time_r(i),1);%right foot in front
    gait_parameters.right.step_length(k) =  markerStruct.RHEE(heel_time_r(i),3)-markerStruct.LHEE(heel_time_r(i),3);%right foot in front
    k=k+1;
end

k=1;
for i= 1:length(heel_time_l)-1
    gait_parameters.left.step_width(k) =  markerStruct.LHEE(heel_time_l(i),1)-markerStruct.RHEE(heel_time_l(i),1);%left foot in front
    gait_parameters.left.step_length(k) =  markerStruct.LHEE(heel_time_l(i),3)-markerStruct.RHEE(heel_time_l(i),3);%left foot in front
    k=k+1;
end
gait_parameters.total.step_length= [gait_parameters.left.step_length,gait_parameters.right.step_length];
gait_parameters.total.step_width= [gait_parameters.left.step_width,gait_parameters.right.step_width];
%% Single Support
k=1;
for i = 2:min(nn_r,nn_l)-1
    if  heel_time_r(i)> heel_time_l(i)
        gait_parameters.left.single_support(k) = (time(heel_time_r(i))-time(toe_time_r(i-1)))/(time(toe_time_l(i))-time(heel_time_l(i)));%percent stance
    else
        gait_parameters.left.single_support(k) = (time(heel_time_r(i+1))-time(toe_time_r(i)))/(time(toe_time_l(i))-time(heel_time_l(i)));%percent stance
    end
    k=k+1;
end
k=1;
for i= 2:min(nn_r,nn_l)-1
    if heel_time_l(i) > heel_time_r(i) 
        gait_parameters.right.single_support(k) = (time(heel_time_l(i))-time(toe_time_l(i-1)))/(time(toe_time_r(i))-time(heel_time_r(i)));%percent stance
    else
        gait_parameters.right.single_support(k) = (time(heel_time_l(i+1)) - time(toe_time_l(i)))/(time(toe_time_r(i)) - time(heel_time_r(i)));%percent stance
    end
    k=k+1;
end

temp_gaitRight = gait_parameters.right.single_support;
temp_gaitLeft = gait_parameters.left.single_support;
gait_parameters.total.single_support = [temp_gaitRight,temp_gaitLeft];

gait_parameters.total.speed = test(50000:end,4)';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

