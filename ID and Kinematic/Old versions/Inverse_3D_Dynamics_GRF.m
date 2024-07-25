%% %%% 
% Code does not include the keeogo output values
%% Static Calibration
% Load OpenSim libs
% clear
% import org.opensim.modeling.*

% % Get the path to a C3D file
% [filename, path] = uigetfile('*.c3d');
% c3dpath = fullfile(path,filename);
% 
% % Construct an opensimC3D object with input c3d path
% % Constructor takes full path to c3d file and an integer for forceplate
% % representation (1 = COP). 
% static = osimC3D(c3dpath,1);
% static.rotateData('x',-90)
% %Get marker and force data
% [static.(name).markers, static.(name).force] = static.getAsStructs();
% a1 = fieldnames(static_markerStruct);
% for i=1:length(a1)-1
%     v = static_markerStruct.(a1{i})(:,1);
%     static_markerStruct.(a1{i})(:,1) = static_markerStruct.(a1{i})(:,3);
%     static_markerStruct.(a1{i})(:,3) = v;
% end
%% EMG
% emg_freq = 2000;
% emg_path = 'C:\Users\MII36\University of Pittsburgh\Keeogo Project - Documents\Pilot Testing\4-27-23\Allison (1)';
% muscle_names = {'lgas','rgas','lham','rham','lquad','rquad'};
% [all_data,MVC,trial_names] = emgAnalysis(emg_path,muscle_names);
% [mvc2] = calculateMVC(MVC.MVC_Knee_ext.rquad);
% % [all_data_norm] = MVCNormalization(all_data,MVC,trial_names,muscle_names);
% [all_data_norm2] = MVCNormalization2(all_data,MVC,trial_names,muscle_names);
% [all_data_norm2] = coactivationratio(all_data_norm2,trial_names,muscle_names);
% time_emg = (1/emg_freq:1/emg_freq:length(all_data.KON_walk.emg.rham)/emg_freq)';
%% Trial
% Get the path to a C3D file
% [filename3, path3] = uigetfile('*.c3d');
% c3dpath = fullfile(path3,filename3);
% c3dExport
function [RKnee, LKnee] = Inverse_3D_Dynamics_GRF(markerStruct, forceStruct, static_markerStruct,anthropometrics,trial_name)
% markerStruct =trial.neutral_control.markers;
% forceStruct =trial.neutral_control.force;
% static_markerStruct = static.neutral_static.markers;
%%
[right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha] = calibrationmatrix_6_27(static_markerStruct,markerStruct);
%[right_globe2ana,left_globe2ana,right_ana2globe,left_ana2globe,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha,RKJC_ana,LKJC_ana] = calibrationmatrix_6_2(static_markerStruct,markerStruct)
% % [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha] = calibrationmatrix2(static_markerStruct,markerStruct);
[foot_right_globe2ana,foot_left_globe2ana,foot_right_seg_omega,foot_left_seg_omega,foot_right_seg_alpha,foot_left_seg_alpha,foot_right_ang,foot_right_omega,foot_right_alpha,foot_left_alpha] = calibrationmatrix_foot_6_27(static_markerStruct,markerStruct);
[right_globe2ana_thigh,left_globe2ana_thigh,right_ana2globe_thigh,left_ana2globe_thigh,right_seg_omega_thigh,left_seg_omega_thigh,right_seg_alpha_thigh,left_seg_alpha_thigh,right_ang_thigh,right_omega_thigh,right_alpha_thigh,left_alpha_thigh,RKJC_ana,LKJC_ana] = calibrationmatrix_thigh(static_markerStruct,markerStruct);
% [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha] = calibrationmatrix_5_25(static_markerStruct,markerStruct);
%%

a = fieldnames(markerStruct);
    % FP filter
    fs = 1000;
    fs_nyq = fs/2;
    win_smooth = 0.25;
    win_size = win_smooth*fs;
    bp_low = 10;
    bp_high = 500;
    fc_low = 25;
    Wn1 = bp_low/fs_nyq;
    Wn2 = bp_high/fs_nyq;
    Wn = fc_low/fs_nyq;
    [bb,ba] = butter(4,Wn,'low');

    % Mocap Filter
    fs_mocap = 100;
    fs_nyq_mo = fs_mocap/2;
    fc_mo_low = 6;
    Wn_mo = fc_mo_low/fs_nyq_mo;
    [cb,ca] = butter(4,Wn_mo,'low');
mass = anthropometrics(1);
height = anthropometrics(2);
gender = anthropometrics(3);

for i = 1:length(a)-1
    markerStruct.(a{i})(isnan(markerStruct.(a{i}))) = 0;
    markerStruct.(a{i}) = filtfilt(cb,ca,markerStruct.(a{i}));
    markerStruct.(a{i}) = markerStruct.(a{i})/1000;
    % m = markerStruct.(a{i})(:,1);
    % markerStruct.(a{i})(:,1) = markerStruct.(a{i})(:,3);
    % markerStruct.(a{i})(:,3) = m;
end

% for i = 1:length(a)-1
%     markerStruct.(a{i}) = filtfilt(cb,ca,markerStruct.(a{i})(:,1:3));
% end
% Find the indices where the condition is met
idx1 = forceStruct.f1(:,2) < 50;
idx2 = forceStruct.f2(:,2) < 50;

% Set the corresponding rows in forceStruct.p1 and forceStruct.m1 to 0
forceStruct.p1(idx1,:) = 0;
forceStruct.m1(idx1,:) = 0;

% Set the corresponding rows in forceStruct.p2 and forceStruct.m2 to 0
forceStruct.p2(idx2,:) = 0;
forceStruct.m2(idx2,:) = 0;

% Update the remaining rows by dividing by 1000
forceStruct.p1(~idx1,:) = forceStruct.p1(~idx1,:) / 1000;
forceStruct.p2(~idx2,:) = forceStruct.p2(~idx2,:) / 1000;
forceStruct.m1(~idx1,:) = forceStruct.m1(~idx1,:) / 1000;
forceStruct.m2(~idx2,:) = forceStruct.m2(~idx2,:) / 1000;

% Filter FP data
forceStruct.m1(:,1) = filtfilt(bb,ba,forceStruct.m1(:,1));
forceStruct.m1(:,2) = filtfilt(bb,ba,forceStruct.m1(:,2));
forceStruct.m1(:,3) = filtfilt(bb,ba,forceStruct.m1(:,3));
forceStruct.m2(:,1) = filtfilt(bb,ba,forceStruct.m2(:,1));
forceStruct.m2(:,2) = filtfilt(bb,ba,forceStruct.m2(:,2));
forceStruct.m2(:,3) = filtfilt(bb,ba,forceStruct.m2(:,3));
forceStruct.f1(:,1:3) = filtfilt(bb,ba,forceStruct.f1(:,1:3));
forceStruct.f2(:,1:3) = filtfilt(bb,ba,forceStruct.f2(:,1:3));
forceStruct.p1(:,1) = filtfilt(bb,ba,forceStruct.p1(:,1));
forceStruct.p1(:,2) = filtfilt(bb,ba,forceStruct.p1(:,2));
forceStruct.p1(:,3) = filtfilt(bb,ba,forceStruct.p1(:,3));
forceStruct.p2(:,1) = filtfilt(bb,ba,forceStruct.p2(:,1));
forceStruct.p2(:,2) = filtfilt(bb,ba,forceStruct.p2(:,2));
forceStruct.p2(:,3) = filtfilt(bb,ba,forceStruct.p2(:,3));
% Filter mocap markers


l_ajc(:,1) = mean([markerStruct.LMedMal(:,1),markerStruct.LLatMal(:,1)],2);
l_ajc(:,2) = mean([markerStruct.LMedMal(:,2),markerStruct.LLatMal(:,2)],2);
l_ajc(:,3) = mean([markerStruct.LMedMal(:,3),markerStruct.LLatMal(:,3)],2);
r_ajc(:,1) = mean([markerStruct.RMedMal(:,1),markerStruct.RLatMal(:,1)],2);
r_ajc(:,2) = mean([markerStruct.RMedMal(:,2),markerStruct.RLatMal(:,2)],2);
r_ajc(:,3) = mean([markerStruct.RMedMal(:,3),markerStruct.RLatMal(:,3)],2);

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

if gender == 1
    Foot.Mass = 0.0137*mass;
    Shank.Mass = 0.0433*mass;
    Foot.MOI(:,1) = Foot.Mass*0.257*(0.258/1.741*height)^2;
    Foot.MOI(:,2) = Foot.Mass*0.245*(0.258/1.741*height)^2;
    Foot.MOI(:,3) = Foot.Mass*0.124*(0.258/1.741*height)^2;
    Shank.MOI(:,1) = Shank.Mass*0.249*(0.434/1.741*height)^2;
    Shank.MOI(:,2) = Shank.Mass*0.103*(0.434/1.741*height)^2;
    Shank.MOI(:,3) = Shank.Mass*0.255*(0.434/1.741*height)^2;
elseif gender == 2
    Foot.Mass = 0.0129*mass;
    Shank.Mass = 0.0481*mass;
    Foot.MOI(:,1) = Foot.Mass*0.299*(0.228/1.735*height)^2;
    Foot.MOI(:,2) = Foot.Mass*0.279*(0.228/1.735*height)^2;
    Foot.MOI(:,3) = Foot.Mass*0.139*(0.228/1.735*height)^2;
    Shank.MOI(:,1) = Shank.Mass*0.267*(0.432/1.735*height)^2;
    Shank.MOI(:,2) = Shank.Mass*0.093*(0.432/1.735*height)^2;
    Shank.MOI(:,3) = Shank.Mass*0.271*(0.432/1.735*height)^2;
end

% Theta.Left.MH(:,1) = unwrap(atan2(markerStruct.LHEE(:,3)-markerStruct.L5Meta(:,3),markerStruct.LHEE(:,2)-markerStruct.L5Meta(:,2))-pi/2);
% Theta.Left.PH(:,1) = unwrap(atan2(markerStruct.L5Meta(:,3)-markerStruct.L2Meta(:,3),markerStruct.L5Meta(:,2)-markerStruct.L2Meta(:,2)));
% Theta.Right.MH(:,1) = unwrap(atan2(markerStruct.RHEE(:,3)-markerStruct.R5Meta(:,3),markerStruct.RHEE(:,2)-markerStruct.R5Meta(:,2)));
% Theta.Right.PH(:,1) = unwrap(atan2(markerStruct.R5Meta(:,3)-markerStruct.R2Meta(:,3),markerStruct.R5Meta(:,2)-markerStruct.R2Meta(:,2)));
% 
% 
% Theta.Left.Shank(:,1) = unwrap(atan2(LKJC_ana(:,3)-l_ajc(:,3),LKJC_ana(:,2)-l_ajc(:,2)));
% Theta.Right.Shank(:,1) = unwrap(atan2(RKJC_ana(:,3)-r_ajc(:,3),RKJC_ana(:,2)-r_ajc(:,2)));
% % Theta.Left.Shank(:,2) = unwrap(atan2(markerStruct.LLatEpi(:,1)-markerStruct.LLatMal(:,1),markerStruct.LLatEpi(:,3)-markerStruct.LLatMal(:,3)));
% Theta.Left.Shank(:,2) = unwrap(atan2(LKJC_ana(:,1)-l_ajc(:,1),LKJC_ana(:,3)-l_ajc(:,3)));
% Theta.Right.Shank(:,2) = unwrap(atan2(RKJC_ana(:,1)-r_ajc(:,1),RKJC_ana(:,3)-r_ajc(:,3)));
% Theta.Left.Shank(:,3) = unwrap(atan2(LKJC_ana(:,2)-l_ajc(:,2),LKJC_ana(:,1)-l_ajc(:,1)));
% Theta.Right.Shank(:,3) = unwrap(atan2(RKJC_ana(:,2)-r_ajc(:,2),RKJC_ana(:,1)-r_ajc(:,1)));
% 
% Theta.Left.Thigh(:,1) = unwrap(atan2(markerStruct.Pelvis1_LeftThigh1_score(:,3)-LKJC_ana(:,3),markerStruct.Pelvis1_LeftThigh1_score(:,2)-LKJC_ana(:,2)));
% Theta.Right.Thigh(:,1) = unwrap(atan2(markerStruct.Pelvis1_RightThigh1_score(:,3)-RKJC_ana(:,3),markerStruct.Pelvis1_RightThigh1_score(:,2)-RKJC_ana(:,2)));
% Theta.Left.Thigh(:,2) = unwrap(atan2(markerStruct.Pelvis1_LeftThigh1_score(:,1)-LKJC_ana(:,1),markerStruct.Pelvis1_LeftThigh1_score(:,3)-LKJC_ana(:,3)));
% Theta.Right.Thigh(:,2) = unwrap(atan2(markerStruct.Pelvis1_RightThigh1_score(:,1)-RKJC_ana(:,1),markerStruct.Pelvis1_RightThigh1_score(:,3)-RKJC_ana(:,3)));
% Theta.Left.Thigh(:,3) = unwrap(atan2(markerStruct.Pelvis1_LeftThigh1_score(:,2)-LKJC_ana(:,2),markerStruct.Pelvis1_LeftThigh1_score(:,1)-LKJC_ana(:,1)));
% Theta.Right.Thigh(:,3) = unwrap(atan2(markerStruct.Pelvis1_RightThigh1_score(:,2)-RKJC_ana(:,2),markerStruct.Pelvis1_RightThigh1_score(:,1)-RKJC_ana(:,1)));

Foot.Left.Theta(:,1) = unwrap(atan2(markerStruct.L5Meta(:,3)-markerStruct.LHEE(:,3),markerStruct.L5Meta(:,2)-markerStruct.LHEE(:,2)));
Foot.Left.Theta(:,2) = unwrap(atan2(markerStruct.L2Meta(:,1)-markerStruct.LHEE(:,1),markerStruct.L2Meta(:,3)-markerStruct.LHEE(:,3)));
Foot.Left.Theta(:,3) = unwrap(atan2(markerStruct.LTOE(:,2)-markerStruct.L5Meta(:,2),markerStruct.LTOE(:,1)-markerStruct.L5Meta(:,1)));
Foot.Right.Theta(:,1) = unwrap(atan2(markerStruct.R5Meta(:,3)-markerStruct.RHEE(:,3),markerStruct.R5Meta(:,2)-markerStruct.RHEE(:,2)));
Foot.Right.Theta(:,2) = unwrap(atan2(markerStruct.R2Meta(:,1)-markerStruct.RHEE(:,1),markerStruct.R2Meta(:,3)-markerStruct.RHEE(:,3)));
Foot.Right.Theta(:,3) = unwrap(atan2(markerStruct.R5Meta(:,2)-markerStruct.RTOE(:,2),markerStruct.R5Meta(:,1)-markerStruct.RTOE(:,1)));

% Theta.Left.Ankle(:,1) = Theta.Left.Shank(:,1)-Theta.Left.MH(:,1) + pi;
% Theta.Left.Ankle(:,2) = unwrap(atan2(markerStruct.L2Meta(:,1)-markerStruct.LHEE(:,1),markerStruct.L2Meta(:,3)-markerStruct.LHEE(:,3)));
% Theta.Left.Ankle(:,3) = unwrap(atan2(markerStruct.LTOE(:,2)-markerStruct.L5Meta(:,2),markerStruct.LTOE(:,1)-markerStruct.L5Meta(:,1)));
% Theta.Right.Ankle(:,1) = Theta.Right.Shank(:,1)-Theta.Right.MH(:,1) + pi;
% Theta.Right.Ankle(:,2) = unwrap(atan2(markerStruct.R2Meta(:,1)-markerStruct.RHEE(:,1),markerStruct.R2Meta(:,3)-markerStruct.RHEE(:,3)));
% Theta.Right.Ankle(:,3) = unwrap(atan2(markerStruct.R5Meta(:,2)-markerStruct.RTOE(:,2),markerStruct.R5Meta(:,1)-markerStruct.RTOE(:,1)));

% LKnee.angle(:,1) = Theta.Left.Thigh(:,1) - Theta.Left.Shank(:,1);
% LKnee.angle(:,2) = Theta.Left.Thigh(:,2) - Theta.Left.Shank(:,2);
% LKnee.angle(:,3) = Theta.Left.Thigh(:,3) - Theta.Left.Shank(:,3);
% 
% RKnee.angle(:,1) = Theta.Right.Thigh(:,1) - Theta.Right.Shank(:,1);
% RKnee.angle(:,2) = Theta.Right.Thigh(:,2) - Theta.Right.Shank(:,2);
% RKnee.angle(:,3) = Theta.Right.Thigh(:,3) - Theta.Right.Shank(:,3);

%% Groud et al. 1983
common_left = cross(left_globe2ana(2,:,:),left_globe2ana_thigh(1,:,:),2)./...
    sqrt(sum(cross(left_globe2ana(2,:,:),left_globe2ana_thigh(1,:,:),2).^2,2));
common_right = cross(right_globe2ana(2,:,:),right_globe2ana_thigh(1,:,:),2)./...
    sqrt(sum(cross(right_globe2ana(2,:,:),right_globe2ana_thigh(1,:,:),2).^2,2));

LKnee.flexion_angle = unwrap(-asin(sum(common_left.*left_globe2ana_thigh(2,:,:),2)));
RKnee.flexion_angle = unwrap(-asin(sum(common_right.*right_globe2ana_thigh(2,:,:),2)));

LKnee.flexion_angle = reshape(LKnee.flexion_angle,length(LKnee.flexion_angle),1);
RKnee.flexion_angle = reshape(RKnee.flexion_angle,length(RKnee.flexion_angle),1);

LKnee.flexion_angle = ScaleTime(LKnee.flexion_angle,1,length(LKnee.flexion_angle),10*length(LKnee.flexion_angle));
RKnee.flexion_angle = ScaleTime(RKnee.flexion_angle,1,length(RKnee.flexion_angle),10*length(RKnee.flexion_angle));

%%

RFoot.COM(isnan(RFoot.COM))= 0;
RFoot.COM = filtfilt(cb,ca,RFoot.COM);
LFoot.COM(isnan(LFoot.COM))= 0;
LFoot.COM = filtfilt(cb,ca,LFoot.COM);
RShank.COM(isnan(RShank.COM))= 0;
RShank.COM = filtfilt(cb,ca,RShank.COM);
LShank.COM(isnan(LShank.COM))= 0;
LShank.COM = filtfilt(cb,ca,LShank.COM);


Velocity.Left.Foot = zeros(size(RFoot.COM));
Velocity.Right.Foot = zeros(size(RFoot.COM));
Velocity.Left.Shank = zeros(size(RFoot.COM));
Velocity.Right.Shank = zeros(size(RFoot.COM));
Acceleration.Left.Foot = zeros(size(RFoot.COM));
Acceleration.Right.Foot = zeros(size(RFoot.COM));
Acceleration.Left.Shank = zeros(size(RFoot.COM));
Acceleration.Right.Shank = zeros(size(RFoot.COM));

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
    % Omega.Left.Ankle(i,1) = (Theta.Left.Ankle(i+1,1) - Theta.Left.Ankle(i-1,1))/.02;
    % Omega.Left.Ankle(i,2) = (Theta.Left.Ankle(i+1,2) - Theta.Left.Ankle(i-1,2))/.02;
    % Omega.Left.Ankle(i,3) = (Theta.Left.Ankle(i+1,3) - Theta.Left.Ankle(i-1,3))/.02;
    % Omega.Right.Ankle(i,1) = (Theta.Right.Ankle(i+1,1) - Theta.Right.Ankle(i-1,1))/.02;
    % Omega.Right.Ankle(i,2) = (Theta.Right.Ankle(i+1,2) - Theta.Right.Ankle(i-1,2))/.02;
    % Omega.Right.Ankle(i,3) = (Theta.Right.Ankle(i+1,3) - Theta.Right.Ankle(i-1,3))/.02;
    % Omega.Left.Shank(i,1) = (Theta.Left.Shank(i+1,1) - Theta.Left.Shank(i-1,1))/.02;
    % Omega.Left.Shank(i,2) = (Theta.Left.Shank(i+1,2) - Theta.Left.Shank(i-1,2))/.02;
    % Omega.Left.Shank(i,3) = (Theta.Left.Shank(i+1,3) - Theta.Left.Shank(i-1,3))/.02;
    % Omega.Right.Shank(i,1) = (Theta.Right.Shank(i+1,1) - Theta.Right.Shank(i-1,1))/.02;
    % Omega.Right.Shank(i,2) = (Theta.Right.Shank(i+1,2) - Theta.Right.Shank(i-1,2))/.02;
    % Omega.Right.Shank(i,3) = (Theta.Right.Shank(i+1,3) - Theta.Right.Shank(i-1,3))/.02;
    % Alpha.Left.Ankle(i,1) = (Theta.Left.Ankle(i+1,1) - 2*Theta.Left.Ankle(i,1) + Theta.Left.Ankle(i-1,1))/(.01^2);
    % Alpha.Left.Ankle(i,2) = (Theta.Left.Ankle(i+1,2) - 2*Theta.Left.Ankle(i,2) + Theta.Left.Ankle(i-1,2))/(.01^2);
    % Alpha.Left.Ankle(i,3) = (Theta.Left.Ankle(i+1,3) - 2*Theta.Left.Ankle(i,3) + Theta.Left.Ankle(i-1,3))/(.01^2);
    % Alpha.Right.Ankle(i,1) = (Theta.Right.Ankle(i+1,1) - 2*Theta.Right.Ankle(i,1) + Theta.Right.Ankle(i-1,1))/(.01^2);
    % Alpha.Right.Ankle(i,2) = (Theta.Right.Ankle(i+1,2) - 2*Theta.Right.Ankle(i,2) + Theta.Right.Ankle(i-1,2))/(.01^2);
    % Alpha.Right.Ankle(i,3) = (Theta.Right.Ankle(i+1,3) - 2*Theta.Right.Ankle(i,3) + Theta.Right.Ankle(i-1,3))/(.01^2);
    % Alpha.Left.Shank(i,1) = (Theta.Left.Shank(i+1,1) - 2*Theta.Left.Shank(i,1) + Theta.Left.Shank(i-1,1))/(.01^2);
    % Alpha.Left.Shank(i,2) = (Theta.Left.Shank(i+1,2) - 2*Theta.Left.Shank(i,2) + Theta.Left.Shank(i-1,2))/(.01^2);
    % Alpha.Left.Shank(i,3) = (Theta.Left.Shank(i+1,3) - 2*Theta.Left.Shank(i,3) + Theta.Left.Shank(i-1,3))/(.01^2);
    % Alpha.Right.Shank(i,1) = (Theta.Right.Shank(i+1,1) - 2*Theta.Right.Shank(i,1) + Theta.Right.Shank(i-1,1))/(.01^2);
    % Alpha.Right.Shank(i,2) = (Theta.Right.Shank(i+1,2) - 2*Theta.Right.Shank(i,2) + Theta.Right.Shank(i-1,2))/(.01^2);
    % Alpha.Right.Shank(i,3) = (Theta.Right.Shank(i+1,3) - 2*Theta.Right.Shank(i,3) + Theta.Right.Shank(i-1,3))/(.01^2);
end
a = fieldnames(Acceleration);
b = fieldnames(Acceleration.Left);
%%
% for i = 1:length(a)
%     for j = 1:length(b)
%         lastValidIndex = find(~isnan(Acceleration.(a{i}).(b{j})),1,'last');
%         trimmedData = Acceleration.(a{i}).(b{j})(1:lastValidIndex);
%         filteredData = filtfilt(cb,ca,trimmedData);
%         Acceleration.(a{i}).(b{j}) = [filteredData; NaN(numel(Acceleration.(a{i}).(b{j})) - numel(filteredData), 3)];
%     end
% end
% for i = 1:length(a)
%     for j = 1:length(b)
%         Acceleration.(a{i}).(b{j})(isnan(Acceleration.(a{i}).(b{j})))= 0;
%         Acceleration.(a{i}).(b{j}) = filtfilt(cb,ca,Acceleration.(a{i}).(b{j}));
%     end
% end
%%
right_globe2ana = ScaleTime_MI(right_globe2ana,1,length(right_globe2ana),10*length(right_globe2ana));
left_globe2ana = ScaleTime_MI(left_globe2ana,1,length(left_globe2ana),10*length(left_globe2ana));
right_globe2ana_thigh = ScaleTime_MI(right_globe2ana_thigh,1,length(right_globe2ana_thigh),10*length(right_globe2ana_thigh));
left_globe2ana_thigh = ScaleTime_MI(left_globe2ana_thigh,1,length(left_globe2ana_thigh),10*length(left_globe2ana_thigh));
RKJC_ana = ScaleTime(RKJC_ana,1,length(RKJC_ana),10*length(RKJC_ana));
LKJC_ana = ScaleTime(LKJC_ana,1,length(LKJC_ana),10*length(LKJC_ana));
a = fieldnames(markerStruct);
for i = 1:length(a)
    markerStruct.(a{i}) = ScaleTime(markerStruct.(a{i}),1,length(markerStruct.(a{i})),10*length(markerStruct.(a{i})));
end
% right_shank_cluster = nanmean([static_markerStruct.RSupLatShank;static_markerStruct.RSupMedShank;...
%     static_markerStruct.RInfLatShank;static_markerStruct.RInfMedShank]);
% left_shank_cluster = nanmean([static_markerStruct.LSupLatShank;static_markerStruct.LSupMedShank;...
%     static_markerStruct.LInfLatShank;static_markerStruct.LInfMedShank]);

LKJC_ana(isnan(LKJC_ana)) = 0;
LKJC_ana = filtfilt(cb,ca,LKJC_ana);
RKJC_ana(isnan(RKJC_ana)) = 0;
RKJC_ana = filtfilt(cb,ca,RKJC_ana);

% right_alpha(isnan(right_alpha)) = 0;
% right_alpha = filtfilt(cb,ca,right_alpha);
% left_alpha(isnan(left_alpha)) = 0;
% left_alpha = filtfilt(cb,ca,left_alpha);
% foot_right_alpha(isnan(foot_right_alpha)) = 0;
% foot_right_alpha = filtfilt(cb,ca,foot_right_alpha);
% foot_left_alpha(isnan(foot_left_alpha)) = 0;
% foot_left_alpha = filtfilt(cb,ca,foot_left_alpha);

right_alpha = ScaleTime(right_alpha,1,length(right_alpha),10*length(right_alpha));
left_alpha = ScaleTime(left_alpha,1,length(left_alpha),10*length(left_alpha));
foot_right_alpha = ScaleTime(foot_right_alpha,1,length(foot_right_alpha),10*length(foot_right_alpha));
foot_left_alpha = ScaleTime(foot_left_alpha,1,length(foot_left_alpha),10*length(foot_left_alpha));
% Alpha.Left.Ankle= ScaleTime(Alpha.Left.Ankle,1,length(Alpha.Left.Ankle),10*length(Alpha.Left.Ankle));
% Alpha.Right.Ankle= ScaleTime(Alpha.Right.Ankle,1,length(Alpha.Right.Ankle),10*length(Alpha.Right.Ankle));
Acceleration.Left.Foot= ScaleTime(Acceleration.Left.Foot,1,length(Acceleration.Left.Foot),10*length(Acceleration.Left.Foot));
Acceleration.Right.Foot= ScaleTime(Acceleration.Right.Foot,1,length(Acceleration.Right.Foot),10*length(Acceleration.Right.Foot));
% Alpha.Left.Shank= ScaleTime(Alpha.Left.Shank,1,length(Alpha.Left.Shank),10*length(Alpha.Left.Shank));
% Alpha.Right.Shank= ScaleTime(Alpha.Right.Shank,1,length(Alpha.Right.Shank),10*length(Alpha.Right.Shank));
Acceleration.Left.Shank= ScaleTime(Acceleration.Left.Shank,1,length(Acceleration.Left.Shank),10*length(Acceleration.Left.Shank));
Acceleration.Right.Shank= ScaleTime(Acceleration.Right.Shank,1,length(Acceleration.Right.Shank),10*length(Acceleration.Right.Shank));
LFoot.COM = ScaleTime(LFoot.COM,1,length(LFoot.COM),10*length(LFoot.COM));
RFoot.COM = ScaleTime(RFoot.COM,1,length(RFoot.COM),10*length(RFoot.COM));
LShank.COM = ScaleTime(LShank.COM,1,length(LShank.COM),10*length(LShank.COM));
RShank.COM = ScaleTime(RShank.COM,1,length(RShank.COM),10*length(RShank.COM));


l_ajc = ScaleTime(l_ajc,1,length(l_ajc),10*length(l_ajc));
r_ajc = ScaleTime(r_ajc,1,length(r_ajc),10*length(r_ajc));
% LAnkle.Force = zeros(3,length(forceStruct.f1))';
% RAnkle.Force = zeros(3,length(forceStruct.f1))';
% LAnkle.Moment = zeros(3,length(forceStruct.f1))';
% RAnkle.Moment = zeros(3,length(forceStruct.f1))';
% 
% time_kin = (1/1000:1/1000:length(forceStruct.f1)/1000)';

% LAnkle.Force(:,1) = Foot.Mass*Acceleration.Left.Foot(:,1) - forceStruct.f1(:,1);
% RAnkle.Force(:,1) = Foot.Mass*Acceleration.Right.Foot(:,1) - forceStruct.f2(:,1);
% if contains(lower(trial_name), 'uphill')
%     LAnkle.Force(:,2) = Foot.Mass*Acceleration.Left.Foot(:,2) - forceStruct.f1(:,2) + Foot.Mass*9.81*cos(8.6*pi/180);
%     RAnkle.Force(:,2) = Foot.Mass*Acceleration.Right.Foot(:,2) - forceStruct.f2(:,2) + Foot.Mass*9.81*cos(8.6*pi/180);
%     LAnkle.Force(:,3) = Foot.Mass*Acceleration.Left.Foot(:,3) - forceStruct.f1(:,3) + Foot.Mass*9.81*sin(8.6*pi/180);
%     RAnkle.Force(:,3) = Foot.Mass*Acceleration.Right.Foot(:,3) - forceStruct.f2(:,3) + Foot.Mass*9.81*sin(8.6*pi/180);
% elseif contains(lower(trial_name), 'down')
%     LAnkle.Force(:,2) = Foot.Mass*Acceleration.Left.Foot(:,2) - forceStruct.f1(:,2) + Foot.Mass*9.81*cos(8.6*pi/180);
%     RAnkle.Force(:,2) = Foot.Mass*Acceleration.Right.Foot(:,2) - forceStruct.f2(:,2) + Foot.Mass*9.81*cos(8.6*pi/180);
%     LAnkle.Force(:,3) = Foot.Mass*Acceleration.Left.Foot(:,3) - forceStruct.f1(:,3) - Foot.Mass*9.81*sin(8.6*pi/180);
%     RAnkle.Force(:,3) = Foot.Mass*Acceleration.Right.Foot(:,3) - forceStruct.f2(:,3) - Foot.Mass*9.81*sin(8.6*pi/180);
% elseif contains(lower(trial_name),'neutral')
%     LAnkle.Force(:,2) = Foot.Mass*Acceleration.Left.Foot(:,2) - forceStruct.f1(:,2) + Foot.Mass*9.81;
%     RAnkle.Force(:,2) = Foot.Mass*Acceleration.Right.Foot(:,2) - forceStruct.f2(:,2) + Foot.Mass*9.81;
%     LAnkle.Force(:,3) = Foot.Mass*Acceleration.Left.Foot(:,3) - forceStruct.f1(:,3);
%     RAnkle.Force(:,3) = Foot.Mass*Acceleration.Right.Foot(:,3) - forceStruct.f2(:,3);
% end

% ind_1 = forceStruct.f1(:,2) < 50;
% ind_2 = forceStruct.f2(:,2) < 50;
% LAnkle.Force(ind_1,:) = 0;
% RAnkle.Force(ind_2,:) = 0;
% Extension
% LAnkle.Moment(:,1) = (Foot.MOI(1))*foot_left_alpha(:,1)-...
%     (forceStruct.f1(:,3) .* (forceStruct.p1(:,2)-LFoot.COM(:,2))) +...
%     (forceStruct.f1(:,2) .* (forceStruct.p1(:,3) - LFoot.COM(:,3))) +...
%     (LAnkle.Force(:,2) .* (l_ajc(:,3)-LFoot.COM(:,3)))-...
%     (LAnkle.Force(:,3) .* (l_ajc(:,2)-LFoot.COM(:,2)));
% %Rotation
% LAnkle.Moment(:,2) = Foot.MOI(2)*foot_left_alpha(:,2)+ ...
%     (forceStruct.f1(:,3) .* (forceStruct.p1(:,1) - LFoot.COM(:,1)))-... 
%     (forceStruct.f1(:,1) .* (forceStruct.p1(:,3) - LFoot.COM(:,3)))+...
%     (LAnkle.Force(:,3) .* (l_ajc(:,1) -LFoot.COM(:,1)))- ...
%     (LAnkle.Force(:,1) .* (l_ajc(:,3)-LFoot.COM(:,3)))-forceStruct.m1(:,2);
% %Abduction
% LAnkle.Moment(:,3) = Foot.MOI(3)*foot_left_alpha(:,3)- ...
%     (forceStruct.f1(:,2) .* (forceStruct.p1(:,1) - LFoot.COM(:,1)))+... 
%     (forceStruct.f1(:,1) .* (forceStruct.p1(:,2) - LFoot.COM(:,2)))-...
%     (LAnkle.Force(:,2) .* (l_ajc(:,1)-LFoot.COM(:,1))) + ...
%     (LAnkle.Force(:,1) .* (l_ajc(:,2)-LFoot.COM(:,2)));
% 
% %Extension
% RAnkle.Moment(:,1) = Foot.MOI(1)*foot_right_alpha(:,1)-...
%     (forceStruct.f2(:,3) .* (forceStruct.p2(:,2)-RFoot.COM(:,2))) +...
%     (forceStruct.f2(:,2) .* (forceStruct.p2(:,3) - RFoot.COM(:,3))) +...
%     (RAnkle.Force(:,2) .* (r_ajc(:,3)-RFoot.COM(:,3)))-...
%     (RAnkle.Force(:,3) .* (r_ajc(:,2)-RFoot.COM(:,2)));
% %Rotation
% RAnkle.Moment(:,2) = Foot.MOI(2)*foot_right_alpha(:,2)+ ...
%     (forceStruct.f2(:,3) .* (forceStruct.p2(:,1) - RFoot.COM(:,1)))-... 
%     (forceStruct.f2(:,1) .* (forceStruct.p2(:,3) - RFoot.COM(:,3)))+...
%     (RAnkle.Force(:,3) .* (r_ajc(:,1)-RFoot.COM(:,1))) - ...
%     (RAnkle.Force(:,1) .* (r_ajc(:,3)-RFoot.COM(:,3)))-forceStruct.m2(:,2);
% %Abduction
% RAnkle.Moment(:,3) = Foot.MOI(3)*foot_right_alpha(:,3)- ...
%     (forceStruct.f2(:,2) .* (forceStruct.p2(:,1) - RFoot.COM(:,1)))+... 
%     (forceStruct.f2(:,1) .* (forceStruct.p2(:,2) - RFoot.COM(:,2)))-...
%     (RAnkle.Force(:,2) .* (r_ajc(:,1)-RFoot.COM(:,1))) + ...
%     (RAnkle.Force(:,1) .* (r_ajc(:,2)-RFoot.COM(:,2)));
% % end

LKnee.Force(:,1) = - forceStruct.f1(:,1);
RKnee.Force(:,1) = - forceStruct.f2(:,1);
LKnee.Force(:,2) =- forceStruct.f1(:,2);
RKnee.Force(:,2) = - forceStruct.f2(:,2);
LKnee.Force(:,3) =- forceStruct.f1(:,3);
RKnee.Force(:,3) = - forceStruct.f2(:,3);


% ind_1 = forceStruct.f1(:,2) < 50;
% ind_2 = forceStruct.f2(:,2) < 50;
% LKnee.Force(ind_1,:) = 0;
% RKnee.Force(ind_2,:) = 0;

% Extension
LKnee.Moment(:,1) = -(forceStruct.f1(:,3) .* (LKJC_ana(:,2))-forceStruct.p1(:,2)) +...
    (forceStruct.f1(:,2) .* (LKJC_ana(:,3)-forceStruct.p1(:,3)));
LKnee.Moment(:,2) = (forceStruct.f1(:,3) .* (LKJC_ana(:,1)-forceStruct.p1(:,1)))-... 
    (forceStruct.f1(:,1) .* (LKJC_ana(:,3)-forceStruct.p1(:,3)))-forceStruct.m1(:,2);
LKnee.Moment(:,3) = -(forceStruct.f1(:,2) .* (LKJC_ana(:,1)-forceStruct.p1(:,1)))+... 
    (forceStruct.f1(:,1) .* (LKJC_ana(:,2)-forceStruct.p1(:,2)));

RKnee.Moment(:,1) = -(forceStruct.f2(:,3) .* (RKJC_ana(:,2)-forceStruct.p2(:,2))) +...
    (forceStruct.f2(:,2) .* (RKJC_ana(:,3)-forceStruct.p2(:,3)));
RKnee.Moment(:,2) = (forceStruct.f2(:,3) .* (RKJC_ana(:,1)-forceStruct.p2(:,1)))-... 
    (forceStruct.f2(:,1) .* (RKJC_ana(:,3)-forceStruct.p2(:,3)))-forceStruct.m2(:,2);
RKnee.Moment(:,3) = -(forceStruct.f2(:,2) .* (RKJC_ana(:,1)-forceStruct.p2(:,1)))+... 
    (forceStruct.f2(:,1) .* (RKJC_ana(:,2)-forceStruct.p2(:,2)));
% end
%% Converting to Segment Coordinate System
% [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha] = calibrationmatrix(static_markerStruct,markerStruct);
    RKnee.right_seg_kneeforce = zeros(length(right_globe2ana),3);
    LKnee.left_seg_kneeforce = zeros(length(right_globe2ana),3);
    % RAnkle.right_seg_anklemom = zeros(length(right_globe2ana),3);
    % LAnkle.left_seg_anklemom = zeros(length(right_globe2ana),3);
    % RAnkle.right_seg_ankleforce = zeros(length(right_globe2ana),3);
    % LAnkle.left_seg_ankleforce = zeros(length(right_globe2ana),3);
    RKnee.right_seg_kneemoment = zeros(length(right_globe2ana),3);
    LKnee.left_seg_kneemoment = zeros(length(right_globe2ana),3);
for i = 1:length(right_globe2ana)
    RKnee.right_seg_kneeforce(i,:) = (right_globe2ana(:,:,i)*(RKnee.Force(i,:)'))';
    LKnee.left_seg_kneeforce(i,:) = (left_globe2ana(:,:,i)*(LKnee.Force(i,:)'))';
    % RAnkle.right_seg_anklemom(i,:) = (right_globe2ana(:,:,i)*(RAnkle.Moment(i,:)'))';
    % LAnkle.left_seg_anklemom(i,:) = (left_globe2ana(:,:,i)*(LAnkle.Moment(i,:)'))';
    % RAnkle.right_seg_ankleforce(i,:) = (right_globe2ana(:,:,i)*(RAnkle.Force(i,:)'))';
    % LAnkle.left_seg_ankleforce(i,:) = (left_globe2ana(:,:,i)*(LAnkle.Force(i,:)'))';
    RKnee.right_seg_kneemoment(i,:) = (right_globe2ana(:,:,i)*(RKnee.Moment(i,:)'))';
    LKnee.left_seg_kneemoment(i,:) = (left_globe2ana(:,:,i)*(LKnee.Moment(i,:)'))';

    RKnee.right_thigh_seg_kneeforce(i,:) = (right_globe2ana_thigh(:,:,i)*(RKnee.Force(i,:)'))';
    LKnee.left_thigh_seg_kneeforce(i,:) = (left_globe2ana_thigh(:,:,i)*(LKnee.Force(i,:)'))';
    RKnee.right_thigh_seg_kneemoment(i,:) = (right_globe2ana_thigh(:,:,i)*(RKnee.Moment(i,:)'))';
    LKnee.left_thigh_seg_kneemoment(i,:) = (left_globe2ana_thigh(:,:,i)*(LKnee.Moment(i,:)'))';
end

% right_seg_kneemoment2 = zeros(length(right_globe2ana),3);
% for i = 1:length(right_globe2ana)
%     if forceStruct.f2(i,2)<25
%         right_seg_kneemoment2(i,1) = 0;
%         right_seg_kneemoment2(i,2) = 0;
%         right_seg_kneemoment2(i,3) = 0;
%     else
%         right_seg_kneemoment2(i,1) = Shank.MOI(1)*right_seg_alpha(i,1) + (Shank.MOI(3)-Shank.MOI(2))*right_seg_omega(i,2)*right_seg_omega(i,3)+...
%         right_seg_anklemom(i,1)+right_seg_kneeforce(i,3)*sqrt(sum((mean([markerStruct.RMedEpi(i,:);markerStruct.RLatEpi(i,:)])-RShank.COM(i,:)).^2))+...
%         (-right_seg_ankleforce(i,3))*sqrt(sum((mean([markerStruct.RMedMal(i,:);markerStruct.RLatMal(i,:)] )-RShank.COM(i,:)).^2));
%     % right_seg_kneemoment2(i,1) = Shank.MOI(1)*right_seg_alpha(i,1) + (Shank.MOI(3)-Shank.MOI(2))*right_seg_omega(i,2)*right_seg_omega(i,3)+...
%     %     right_seg_anklemom(i,1)+right_seg_kneeforce(i,3)*sqrt((mean([markerStruct.RMedEpi(i,1),markerStruct.RLatEpi(i,1)]-RShank.COM(i,1))^2)+...
%     %     (mean([markerStruct.RMedEpi(i,2),markerStruct.RLatEpi(i,2)]-RShank.COM(i,2))^2)+...
%     %     (mean([markerStruct.RMedEpi(i,3),markerStruct.RLatEpi(i,3)]-RShank.COM(i,3))^2))+...
%     %     right_seg_ankleforce(i,3)*sqrt((mean([markerStruct.RMedMal(i,1),markerStruct.RLatMal(i,1)]-RShank.COM(i,1))^2)+...
%     %     (mean([markerStruct.RMedMal(i,2),markerStruct.RLatMal(i,2)]-RShank.COM(i,2))^2)+...
%     %     (mean([markerStruct.RMedMal(i,3),markerStruct.RLatMal(i,3)]-RShank.COM(i,3))^2));
%     right_seg_kneemoment2(i,2) = Shank.MOI(2)*right_seg_alpha(i,2) + (Shank.MOI(1)-Shank.MOI(3))*right_seg_omega(i,1)*right_seg_omega(i,3)+...
%         right_seg_anklemom(i,2);
%     right_seg_kneemoment2(i,3) = Shank.MOI(3)*right_seg_alpha(i,3) + (Shank.MOI(2)-Shank.MOI(1))*right_seg_omega(i,1)*right_seg_omega(i,2)+...
%         right_seg_anklemom(i,3)-right_seg_kneeforce(i,1)*sqrt((mean([markerStruct.RMedEpi(i,2),markerStruct.RLatEpi(i,1)]-RShank.COM(i,1))^2)+...
%         (mean([markerStruct.RMedEpi(i,2),markerStruct.RLatEpi(i,2)]-RShank.COM(i,2))^2)+...
%         (mean([markerStruct.RMedEpi(i,3),markerStruct.RLatEpi(i,3)]-RShank.COM(i,3))^2))+...
%         right_seg_ankleforce(i,1)*sqrt((mean([markerStruct.RMedMal(i,1),markerStruct.RLatMal(i,1)]-RShank.COM(i,1))^2)+...
%         (mean([markerStruct.RMedMal(i,2),markerStruct.RLatMal(i,2)]-RShank.COM(i,2))^2)+...
%         (mean([markerStruct.RMedMal(i,3),markerStruct.RLatMal(i,3)]-RShank.COM(i,3))^2));
%     end
% end
% Step_Cycle_MI
end