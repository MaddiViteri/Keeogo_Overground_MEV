%% %%% 
% Test code
% markerStruct = trial.neutral_control.markers;
% forceStruct = trial.neutral_control.force;
% static_markerStruct = static.neutral_static.markers;

% Code does not include the keeogo output values
function [Knee, Ankle, markerStruct, forceStruct] = Inverse_3D_Dynamics_COM_7_18(markerStruct, forceStruct, static_markerStruct,anthropometrics,trial_name)

%% Determine Segment Orientation
% Shank
[right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,...
    right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,...
    left_alpha] = calibrationmatrix_6_27(static_markerStruct,markerStruct);

% Foot
[foot_right_globe2ana,foot_left_globe2ana,foot_right_seg_omega,...
    foot_left_seg_omega,foot_right_seg_alpha,foot_left_seg_alpha,...
    foot_right_ang,foot_right_omega,foot_right_alpha,foot_left_alpha] ...
    = calibrationmatrix_foot_6_27(static_markerStruct,markerStruct);

% Thigh
[right_globe2ana_thigh,left_globe2ana_thigh,right_ana2globe_thigh,...
    left_ana2globe_thigh,right_seg_omega_thigh,left_seg_omega_thigh,...
    right_seg_alpha_thigh,left_seg_alpha_thigh,right_ang_thigh,...
    right_omega_thigh,right_alpha_thigh,left_alpha_thigh,RKJC_ana,LKJC_ana]...
    = calibrationmatrix_thigh(static_markerStruct,markerStruct);

%%

    % Mocap Filter
    fs_mocap = 100;
    fs_nyq_mo = fs_mocap/2;
    fc_mo_low = 6;
    Wn_mo = fc_mo_low/fs_nyq_mo;
    [cb,ca] = butter(4,Wn_mo,'low');

mass = anthropometrics(1);
height = anthropometrics(2);
gender = anthropometrics(3);

% Filter mocap markers
a = fieldnames(markerStruct);
for i = 1:length(a)-1

    markerStruct.(a{i})(isnan(markerStruct.(a{i}))) = 0;
    markerStruct.(a{i}) = filtfilt(cb,ca,markerStruct.(a{i}));
    markerStruct.(a{i}) = markerStruct.(a{i})/1000;

end

%Find Ankle Joint Centers
for i = 1:3
    l_ajc(:,i) = mean([markerStruct.LMedMal(:,i), markerStruct.LLatMal(:,i)], 2);
    r_ajc(:,i) = mean([markerStruct.RMedMal(:,i), markerStruct.RLatMal(:,i)], 2);
end

% Find Segment Center of Mass Using 
LFoot.COM(:,1) = (markerStruct.LTOE(:,1) + markerStruct.L5Meta(:,1))*0.5;
LFoot.COM(:,2) = (markerStruct.LLatMal(:,2) + markerStruct.L2Meta(:,2))*0.5;
LFoot.COM(:,3) = (markerStruct.LLatMal(:,3) + markerStruct.L2Meta(:,3))*0.5;
RFoot.COM(:,1) = (markerStruct.RTOE(:,1) + markerStruct.R5Meta(:,1))*0.5;
RFoot.COM(:,2) = (markerStruct.RLatMal(:,2) + markerStruct.R2Meta(:,2))*0.5;
RFoot.COM(:,3) = (markerStruct.RLatMal(:,3) + markerStruct.R2Meta(:,3))*0.5;
LShank.COM(:,1) = (0.433*markerStruct.LMedMal(:,1) + 0.567*LKJC_ana(:,1));
LShank.COM(:,2) = (0.433*markerStruct.LMedMal(:,2) + 0.567*LKJC_ana(:,2));
LShank.COM(:,3) = (0.433*markerStruct.LMedMal(:,3) + 0.567*LKJC_ana(:,3));
RShank.COM(:,1) = (0.433*markerStruct.RMedMal(:,1) + 0.567*RKJC_ana(:,1));
RShank.COM(:,2) = (0.433*markerStruct.RMedMal(:,2) + 0.567*RKJC_ana(:,2));
RShank.COM(:,3) = (0.433*markerStruct.RMedMal(:,3) + 0.567*RKJC_ana(:,3));

%% Grood and Suntay 1983
% Angle Calculation
common_left = cross(left_globe2ana(2,:,:),left_globe2ana_thigh(1,:,:),2)./...
    sqrt(sum(cross(left_globe2ana(2,:,:),left_globe2ana_thigh(1,:,:),2).^2,2));
common_right = cross(right_globe2ana(2,:,:),right_globe2ana_thigh(1,:,:),2)./...
    sqrt(sum(cross(right_globe2ana(2,:,:),right_globe2ana_thigh(1,:,:),2).^2,2));

Knee.Left.flexion_angle = unwrap(-asin(sum(common_left.*left_globe2ana_thigh(2,:,:),2)));
Knee.Right.flexion_angle = unwrap(-asin(sum(common_right.*right_globe2ana_thigh(2,:,:),2)));

%Reshape the table
Knee.Left.flexion_angle = reshape(Knee.Left.flexion_angle,length(Knee.Left.flexion_angle),1);
Knee.Right.flexion_angle = reshape(Knee.Right.flexion_angle,length(Knee.Right.flexion_angle),1);

%Upsample data to match force plate
Knee.Left.flexion_angle = ScaleTime(Knee.Left.flexion_angle,1,length(Knee.Left.flexion_angle),10*length(Knee.Left.flexion_angle));
Knee.Right.flexion_angle = ScaleTime(Knee.Right.flexion_angle,1,length(Knee.Right.flexion_angle),10*length(Knee.Right.flexion_angle));

% Ankle
common_left_ank = cross(foot_left_globe2ana(2,:,:),left_globe2ana(1,:,:),2)./...
    sqrt(sum(cross(foot_left_globe2ana(2,:,:),left_globe2ana(1,:,:),2).^2,2));
common_right_ank = cross(foot_right_globe2ana(2,:,:),right_globe2ana(1,:,:),2)./...
    sqrt(sum(cross(foot_right_globe2ana(2,:,:),right_globe2ana(1,:,:),2).^2,2));

Ankle.Left.flexion_angle = unwrap(-asin(sum(common_left_ank.*left_globe2ana(2,:,:),2)));
Ankle.Right.flexion_angle = unwrap(-asin(sum(common_right_ank.*right_globe2ana(2,:,:),2)));

%Reshape the table
Ankle.Left.flexion_angle = reshape(Ankle.Left.flexion_angle,length(Ankle.Left.flexion_angle),1);
Ankle.Right.flexion_angle = reshape(Ankle.Right.flexion_angle,length(Ankle.Right.flexion_angle),1);

%Upsample data to match force plate
Ankle.Left.flexion_angle = ScaleTime(Ankle.Left.flexion_angle,1,length(Ankle.Left.flexion_angle),10*length(Ankle.Left.flexion_angle));
Ankle.Right.flexion_angle = ScaleTime(Ankle.Right.flexion_angle,1,length(Ankle.Right.flexion_angle),10*length(Ankle.Right.flexion_angle));

%%
% Filter data
RFoot.COM(isnan(RFoot.COM))= 0;
RFoot.COM = filtfilt(cb,ca,RFoot.COM);
LFoot.COM(isnan(LFoot.COM))= 0;
LFoot.COM = filtfilt(cb,ca,LFoot.COM);

RShank.COM(isnan(RShank.COM))= 0;
RShank.COM = filtfilt(cb,ca,RShank.COM);
LShank.COM(isnan(LShank.COM))= 0;
LShank.COM = filtfilt(cb,ca,LShank.COM); 

LKJC_ana(isnan(LKJC_ana)) = 0;
LKJC_ana = filtfilt(cb,ca,LKJC_ana);
RKJC_ana(isnan(RKJC_ana)) = 0;
RKJC_ana = filtfilt(cb,ca,RKJC_ana);

% Initialize Velocity and Acceleration
Velocity.Left.Foot = zeros(size(RFoot.COM));
Velocity.Right.Foot = zeros(size(RFoot.COM));
Velocity.Left.Shank = zeros(size(RFoot.COM));
Velocity.Right.Shank = zeros(size(RFoot.COM));
Acceleration.Left.Foot = zeros(size(RFoot.COM));
Acceleration.Right.Foot = zeros(size(RFoot.COM));
Acceleration.Left.Shank = zeros(size(RFoot.COM));
Acceleration.Right.Shank = zeros(size(RFoot.COM));

%Calculate Acceleration and Velocity with Position Data
for i = 2:length(RFoot.COM)-1
    Velocity.Left.Foot(i,1) = (LFoot.COM(i+1,1) - LFoot.COM(i-1,1))/.02;
    Velocity.Left.Foot(i,2) = (LFoot.COM(i+1,2) - LFoot.COM(i-1,2))/.02; 
    Velocity.Left.Foot(i,3) = (LFoot.COM(i+1,3) - LFoot.COM(i-1,3))/.02; 
    Velocity.Right.Foot(i,1) = (RFoot.COM(i+1,1) - RFoot.COM(i-1,1))/.02; 
    Velocity.Right.Foot(i,2) = (RFoot.COM(i+1,2) - RFoot.COM(i-1,2))/.02; 
    Velocity.Right.Foot(i,3) = (RFoot.COM(i+1,3) - RFoot.COM(i-1,3))/.02; 
    Velocity.Left.Shank(i,1) = (LShank.COM(i+1,1) - LShank.COM(i-1,1))/.02; 
    Velocity.Left.Shank(i,2) = (LShank.COM(i+1,2) - LShank.COM(i-1,2))/.02; 
    Velocity.Left.Shank(i,3) = (LShank.COM(i+1,3) - LShank.COM(i-1,3))/.02; 
    Velocity.Right.Shank(i,1) = (RShank.COM(i+1,1) - RShank.COM(i-1,1))/.02; 
    Velocity.Right.Shank(i,2) = (RShank.COM(i+1,2) - RShank.COM(i-1,2))/.02; 
    Velocity.Right.Shank(i,3) = (RShank.COM(i+1,3) - RShank.COM(i-1,3))/.02; 
    Acceleration.Left.Foot(i,1) = (LFoot.COM(i+1,1) - 2*LFoot.COM(i,1) + LFoot.COM(i-1,1))/(.01^2);
    Acceleration.Left.Foot(i,2) = (LFoot.COM(i+1,2) - 2*LFoot.COM(i,2) + LFoot.COM(i-1,2))/(.01^2); 
    Acceleration.Left.Foot(i,3) = (LFoot.COM(i+1,3) - 2*LFoot.COM(i,3) + LFoot.COM(i-1,3))/(.01^2); 
    Acceleration.Right.Foot(i,1) = (RFoot.COM(i+1,1) - 2*RFoot.COM(i,1) + RFoot.COM(i-1,1))/(.01^2); 
    Acceleration.Right.Foot(i,2) = (RFoot.COM(i+1,2) - 2*RFoot.COM(i,2) + RFoot.COM(i-1,2))/(.01^2); 
    Acceleration.Right.Foot(i,3) = (RFoot.COM(i+1,3) - 2*RFoot.COM(i,3) + RFoot.COM(i-1,3))/(.01^2); 
    Acceleration.Left.Shank(i,1) = (LShank.COM(i+1,1) - 2*LShank.COM(i,1) + LShank.COM(i-1,1))/(.01^2); 
    Acceleration.Left.Shank(i,2) = (LShank.COM(i+1,2) - 2*LShank.COM(i,2) + LShank.COM(i-1,2))/(.01^2); 
    Acceleration.Left.Shank(i,3) = (LShank.COM(i+1,3) - 2*LShank.COM(i,3) + LShank.COM(i-1,3))/(.01^2); 
    Acceleration.Right.Shank(i,1) = (RShank.COM(i+1,1) - 2*RShank.COM(i,1) + RShank.COM(i-1,1))/(.01^2);
    Acceleration.Right.Shank(i,2) = (RShank.COM(i+1,2) - 2*RShank.COM(i,2) + RShank.COM(i-1,2))/(.01^2);
    Acceleration.Right.Shank(i,3) = (RShank.COM(i+1,3) - 2*RShank.COM(i,3) + RShank.COM(i-1,3))/(.01^2);
end

%%
% Upsample data to match with force plate data
right_globe2ana = ScaleTime_MI(right_globe2ana,1,length(right_globe2ana),10*length(right_globe2ana));
left_globe2ana = ScaleTime_MI(left_globe2ana,1,length(left_globe2ana),10*length(left_globe2ana));
right_globe2ana_thigh = ScaleTime_MI(right_globe2ana_thigh,1,length(right_globe2ana_thigh),10*length(right_globe2ana_thigh));
left_globe2ana_thigh = ScaleTime_MI(left_globe2ana_thigh,1,length(left_globe2ana_thigh),10*length(left_globe2ana_thigh));
RKJC_ana = ScaleTime(RKJC_ana,1,length(RKJC_ana),10*length(RKJC_ana));
LKJC_ana = ScaleTime(LKJC_ana,1,length(LKJC_ana),10*length(LKJC_ana));

% Upsample marker data
a = fieldnames(markerStruct);
for i = 1:length(a)
    markerStruct.(a{i}) = ScaleTime(markerStruct.(a{i}),1,length(markerStruct.(a{i})),10*length(markerStruct.(a{i})));
end

% Upsample Shank Angular Acceleration
right_alpha = ScaleTime(right_alpha,1,length(right_alpha),10*length(right_alpha));
left_alpha = ScaleTime(left_alpha,1,length(left_alpha),10*length(left_alpha));

% Upsample Foot Angular Acceleration
foot_right_alpha = ScaleTime(foot_right_alpha,1,length(foot_right_alpha),10*length(foot_right_alpha));
foot_left_alpha = ScaleTime(foot_left_alpha,1,length(foot_left_alpha),10*length(foot_left_alpha));

% Upsample Foot Acceleration
Acceleration.Left.Foot= ScaleTime(Acceleration.Left.Foot,1,length(Acceleration.Left.Foot),10*length(Acceleration.Left.Foot));
Acceleration.Right.Foot= ScaleTime(Acceleration.Right.Foot,1,length(Acceleration.Right.Foot),10*length(Acceleration.Right.Foot));

% Upsample Shank Acceleration
Acceleration.Left.Shank= ScaleTime(Acceleration.Left.Shank,1,length(Acceleration.Left.Shank),10*length(Acceleration.Left.Shank));
Acceleration.Right.Shank= ScaleTime(Acceleration.Right.Shank,1,length(Acceleration.Right.Shank),10*length(Acceleration.Right.Shank));

% Upsample Segment Position
LFoot.COM = ScaleTime(LFoot.COM,1,length(LFoot.COM),10*length(LFoot.COM));
RFoot.COM = ScaleTime(RFoot.COM,1,length(RFoot.COM),10*length(RFoot.COM));
LShank.COM = ScaleTime(LShank.COM,1,length(LShank.COM),10*length(LShank.COM));
RShank.COM = ScaleTime(RShank.COM,1,length(RShank.COM),10*length(RShank.COM));

% Upsample Ankle Joint Position
l_ajc = ScaleTime(l_ajc,1,length(l_ajc),10*length(l_ajc));
r_ajc = ScaleTime(r_ajc,1,length(r_ajc),10*length(r_ajc));


end