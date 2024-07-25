%% Static Calibration Foot
% Calculates the euler angles and rotation matrix for the foot
% Load OpenSim libs
function [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha] = calibrationmatrix_foot_6_27(static_markerStruct,markerStruct)
%% Marker to anatomical  
% Set filters
fs_mocap = 100;
    fs_nyq_mo = fs_mocap/2;
    fc_mo_low = 6;
    Wn_mo = fc_mo_low/fs_nyq_mo;
    [cb,ca] = butter(4,Wn_mo,'low');

% Determine center of mass from static trial 
left.static_COM(:,1) = (mean(static_markerStruct.LTOE(:,1),"omitnan") + mean(static_markerStruct.L5Meta(:,1),"omitnan"))*0.5;
left.static_COM(:,2) = (mean(static_markerStruct.LLatMal(:,2),"omitnan") + mean(static_markerStruct.L2Meta(:,2),"omitnan"))*0.5;
left.static_COM(:,3) = (mean(static_markerStruct.LLatMal(:,3),"omitnan") + mean(static_markerStruct.L2Meta(:,3),"omitnan"))*0.5;
right.static_COM(:,1) = (mean(static_markerStruct.RTOE(:,1),"omitnan") + mean(static_markerStruct.R5Meta(:,1),"omitnan"))*0.5;
right.static_COM(:,2) = (mean(static_markerStruct.RLatMal(:,2),"omitnan") + mean(static_markerStruct.R2Meta(:,2),"omitnan"))*0.5;
right.static_COM(:,3) = (mean(static_markerStruct.RLatMal(:,3),"omitnan") + mean(static_markerStruct.R2Meta(:,3),"omitnan"))*0.5;

%% Anatomical (an) Axes calculation
% Right Foot
right_foot_xan = (mean(static_markerStruct.RTOE,"omitnan") - mean(static_markerStruct.R5Meta,"omitnan"))/...
    norm((mean(static_markerStruct.RTOE,"omitnan") - mean(static_markerStruct.R5Meta,"omitnan")));
right_foot_zan = (mean(static_markerStruct.R2Meta,"omitnan") - mean(static_markerStruct.RHEE,"omitnan"))/...
    norm((mean(static_markerStruct.R2Meta,"omitnan") - mean(static_markerStruct.RHEE,"omitnan")));
right_foot_yan = cross(right_foot_zan,right_foot_xan);
right_foot_xan = cross(right_foot_yan,right_foot_zan);
right_foot_xan = right_foot_xan/norm(right_foot_xan);
right_foot_yan = right_foot_yan/norm(right_foot_yan);
right_foot_zan = right_foot_zan/norm(right_foot_zan);

% Left Foot
left_foot_xan = (mean(static_markerStruct.L5Meta,"omitnan") - mean(static_markerStruct.LTOE,"omitnan"))/...
    norm((mean(static_markerStruct.L5Meta,"omitnan") - mean(static_markerStruct.LTOE,"omitnan")));
left_foot_zan = (mean(static_markerStruct.L2Meta,"omitnan") - mean(static_markerStruct.LHEE,"omitnan"))/...
    norm((mean(static_markerStruct.L2Meta,"omitnan") - mean(static_markerStruct.LHEE,"omitnan")));
left_foot_yan = cross(left_foot_zan,left_foot_xan);
left_foot_xan = cross(left_foot_yan,left_foot_zan);
left_foot_xan = left_foot_xan/norm(left_foot_xan);
left_foot_yan = left_foot_yan/norm(left_foot_yan);
left_foot_zan = left_foot_zan/norm(left_foot_zan);

% Matrices
RLG2A = [right_foot_xan;right_foot_yan;right_foot_zan];
LLG2A = [left_foot_xan;left_foot_yan;left_foot_zan];
RLA2G = RLG2A';
LLA2G = LLG2A';

%% Marker's relative anatomical vector position
%% 2 Meta
% Right
right_foot2meta2_xal = mean(static_markerStruct.R2Meta(:,1),"omitnan") - right.static_COM(1);
right_foot2meta2_yal = mean(static_markerStruct.R2Meta(:,2),"omitnan") - right.static_COM(2);
right_foot2meta2_zal = mean(static_markerStruct.R2Meta(:,3),"omitnan") - right.static_COM(3);
right_foot2meta2_al = [right_foot2meta2_xal;right_foot2meta2_yal;right_foot2meta2_zal];
right_an_2meta = RLG2A*right_foot2meta2_al;

% Left
left_foot2meta2_xal = mean(static_markerStruct.L2Meta(:,1),"omitnan") - left.static_COM(1);
left_foot2meta2_yal = mean(static_markerStruct.L2Meta(:,2),"omitnan") - left.static_COM(2);
left_foot2meta2_zal = mean(static_markerStruct.L2Meta(:,3),"omitnan") - left.static_COM(3);
left_foot2meta2_al = [left_foot2meta2_xal;left_foot2meta2_yal;left_foot2meta2_zal];
left_an_2meta = LLG2A*left_foot2meta2_al;

%% Heel
% Right
right_foot2heel_xal = mean(static_markerStruct.RHEE(:,1),"omitnan") - right.static_COM(1);
right_foot2heel_yal = mean(static_markerStruct.RHEE(:,2),"omitnan") - right.static_COM(2);
right_foot2heel_zal = mean(static_markerStruct.RHEE(:,3),"omitnan") - right.static_COM(3);
right_foot2heel_al = [right_foot2heel_xal;right_foot2heel_yal;right_foot2heel_zal];
right_an_heel = RLG2A*right_foot2heel_al;
% Left
left_foot2heel_xal = mean(static_markerStruct.LHEE(:,1),"omitnan") - left.static_COM(1);
left_foot2heel_yal = mean(static_markerStruct.LHEE(:,2),"omitnan") - left.static_COM(2);
left_foot2heel_zal = mean(static_markerStruct.LHEE(:,3),"omitnan") - left.static_COM(3);
left_foot2heel_al = [left_foot2heel_xal;left_foot2heel_yal;left_foot2heel_zal];
left_an_heel = LLG2A*left_foot2heel_al;

%% 5th Meta
% Right
right_foot2meta5_xal = mean(static_markerStruct.R5Meta(:,1),"omitnan") - right.static_COM(1);
right_foot2meta5_yal = mean(static_markerStruct.R5Meta(:,2),"omitnan") - right.static_COM(2);
right_foot2meta5_zal = mean(static_markerStruct.R5Meta(:,3),"omitnan") - right.static_COM(3);
right_foot2meta5_al = [right_foot2meta5_xal;right_foot2meta5_yal;right_foot2meta5_zal];
right_an_5meta = RLG2A*right_foot2meta5_al;

% Left
left_foot2meta5_xal = mean(static_markerStruct.L5Meta(:,1),"omitnan") - left.static_COM(1);
left_foot2meta5_yal = mean(static_markerStruct.L5Meta(:,2),"omitnan") - left.static_COM(2);
left_foot2meta5_zal = mean(static_markerStruct.L5Meta(:,3),"omitnan") - left.static_COM(3);
left_foot2meta5_al = [left_foot2meta5_xal;left_foot2meta5_yal;left_foot2meta5_zal];
left_an_5meta = LLG2A*left_foot2meta5_al;

%% TOE
% Right
right_foot2toe_xal = mean(static_markerStruct.RTOE(:,1),"omitnan") - right.static_COM(1);
right_foot2toe_yal = mean(static_markerStruct.RTOE(:,2),"omitnan") - right.static_COM(2);
right_foot2toe_zal = mean(static_markerStruct.RTOE(:,3),"omitnan") - right.static_COM(3);
right_foot2toe_al = [right_foot2toe_xal;right_foot2toe_yal;right_foot2toe_zal];
right_an_toe = RLG2A*right_foot2toe_al;
% Left
left_foot2LRfoot_xal = mean(static_markerStruct.LTOE(:,1),"omitnan") - left.static_COM(1);
left_foot2LRfoot_yal = mean(static_markerStruct.LTOE(:,2),"omitnan") - left.static_COM(2);
left_foot2LRfoot_zal = mean(static_markerStruct.LTOE(:,3),"omitnan") - left.static_COM(3);
left_foot2LRfoot_al = [left_foot2LRfoot_xal;left_foot2LRfoot_yal;left_foot2LRfoot_zal];
left_an_toe = LLG2A*left_foot2LRfoot_al;

%% Marker Matrix
% Right
right_local_xm = (right_an_toe - right_an_heel)/norm(right_an_toe - right_an_heel);
right_local_a = (right_an_2meta - right_an_heel)/norm(right_an_2meta - right_an_heel);
right_local_ym = cross(right_local_xm,right_local_a);
right_local_zm = cross(right_local_xm,right_local_ym);
right_local_zm = right_local_zm/norm(right_local_zm);
right_local_ym = right_local_ym/norm(right_local_ym);
right_local_xm = right_local_xm/norm(right_local_xm);
right_LM2A = [right_local_xm,right_local_ym,right_local_zm];

% Left
left_local_xm = (left_an_toe - left_an_heel)/norm(left_an_toe - left_an_heel);
left_local_a = (left_an_2meta - left_an_heel)/norm(left_an_toe - left_an_heel);
left_local_ym = cross(left_local_xm,left_local_a);
left_local_zm = cross(left_local_xm,left_local_ym);
left_local_zm = left_local_zm/norm(left_local_zm);
left_local_ym = left_local_ym/norm(left_local_ym);
left_local_xm = left_local_xm/norm(left_local_xm);
left_LM2A = [left_local_xm,left_local_ym,left_local_zm];

%% Global to Marker
% Create rotation matrix to get from global positioning to marker position
% as mediary step

% Initialize Matrices
right_globe2mark = zeros(3,3,length(markerStruct.RLatMal));
left_globe2mark = zeros(3,3,length(markerStruct.RLatMal));

% Right
right_globe_xm = (markerStruct.RTOE - markerStruct.RHEE)./...
    sqrt(sum((markerStruct.RTOE - markerStruct.RHEE).^2,2));
right_globe_am = (markerStruct.R2Meta - markerStruct.RHEE)./...
    sqrt(sum((markerStruct.R2Meta - markerStruct.RHEE).^2,2));
right_globe_ym = cross(right_globe_xm,right_globe_am,2);
right_globe_ym = right_globe_ym./sqrt(sum(right_globe_ym.^2,2));
right_globe_zm = cross(right_globe_xm, right_globe_ym,2);
right_globe_zm = right_globe_zm./sqrt(sum(right_globe_zm.^2,2));

% Left
left_globe_xm = (markerStruct.LTOE - markerStruct.LHEE)./...
    sqrt(sum((markerStruct.LTOE - markerStruct.LHEE).^2,2));
left_globe_am = (markerStruct.L2Meta - markerStruct.LHEE)./...
    sqrt(sum((markerStruct.L2Meta - markerStruct.LHEE).^2,2));
left_globe_ym = cross(left_globe_xm,left_globe_am,2);
left_globe_ym = left_globe_ym./sqrt(sum(left_globe_ym.^2,2));
left_globe_zm = cross(left_globe_xm,left_globe_ym,2);
left_globe_zm = left_globe_zm./sqrt(sum(left_globe_zm.^2,2));

for i = 1:length(left_globe_zm)
    right_globe2mark(:,:,i) =  [right_globe_xm(i,:);right_globe_ym(i,:);right_globe_zm(i,:)];
    left_globe2mark(:,:,i) =  [left_globe_xm(i,:);left_globe_ym(i,:);left_globe_zm(i,:)];
end
%% Global to Anatomical
% Get from the Global from the mediary step and matrix multiplication with
% the static calibration that generated the relationship between markers
% and anatomical

for i = 1:length(right_globe2mark)
    right_globe2ana(:,:,i) = right_LM2A*right_globe2mark(:,:,i);
    left_globe2ana(:,:,i) = left_LM2A*left_globe2mark(:,:,i);
end

%% Angle Calculation
% General premise of calculating the joint/segment angles
% syms x y z
% mat1 = [1 0 0;0 cos(x) sin(x); 0 -sin(x) cos(x)];
% mat2 = [cos(y) 0 -sin(y);0 1 0;-sin(y) 0 cos(y)];
% mat3 = [cos(z) sin(z) 0;-sin(z) cos(z) 0;0 0 1];
% rot = mat3*mat2*mat1;

% Calculate angles
right_ang_alpha = unwrap(atan2(-right_globe2ana(3,2,:),right_globe2ana(3,3,:)));
right_ang_beta = unwrap(atan2(right_globe2ana(3,1,:),sqrt((right_globe2ana(1,1,:)).^2+(right_globe2ana(2,1,:)).^2)));
right_ang_gamma = unwrap(atan2(-right_globe2ana(2,1,:),right_globe2ana(1,1,:)));
right_ang = [right_ang_alpha(:),right_ang_beta(:),right_ang_gamma(:)];
left_ang_alpha = unwrap(atan2(-left_globe2ana(3,2,:),left_globe2ana(3,3,:)));
left_ang_beta = unwrap(atan2(left_globe2ana(3,1,:),sqrt((left_globe2ana(1,1,:)).^2+(left_globe2ana(2,1,:)).^2)));
left_ang_gamma = unwrap(atan2(-left_globe2ana(2,1,:),left_globe2ana(1,1,:)));
left_ang = [left_ang_alpha(:),left_ang_beta(:),left_ang_gamma(:)];

% Filter angular position
right_ang(isnan(right_ang)) = 0;
right_ang = filtfilt(cb,ca,right_ang);
left_ang(isnan(left_ang)) = 0;
left_ang = filtfilt(cb,ca,left_ang);

% Initialize Variables for omega and alpha
right_omega = zeros(length(right_ang_alpha),3);
left_omega = zeros(length(left_ang_alpha),3);
right_alpha = zeros(length(right_ang_alpha),3);
left_alpha = zeros(length(left_ang_alpha),3);

% Calculate angular velocity and acceleration 
for i = 2:length(right_ang_alpha)-1
    right_omega(i,1) = (right_ang(i+1,1)-right_ang(i-1,1))/.02;
    right_omega(i,2) = (right_ang(i+1,2)-right_ang(i-1,2))/.02;
    right_omega(i,3) = (right_ang(i+1,3)-right_ang(i-1,3))/.02;
    left_omega(i,1) = (left_ang(i+1,1)-left_ang(i-1,1))/.02;
    left_omega(i,2) = (left_ang(i+1,2)-left_ang(i-1,2))/.02;
    left_omega(i,3) = (left_ang(i+1,3)-left_ang(i-1,3))/.02;
    right_alpha(i,1) = (right_ang(i+1,1) - 2*right_ang(i,1) + right_ang(i-1,1))/(.01^2);
    right_alpha(i,2) = (right_ang(i+1,2)-2*right_ang(i,2) +right_ang(i-1,2))/(.01^2);
    right_alpha(i,3) = (right_ang(i+1,3)-2*right_ang(i,3)+right_ang(i-1,3))/(.01^2);
    left_alpha(i,1) = (left_ang(i+1,1)-2*left_ang(i,1)+left_ang(i-1,1))/(.01^2);
    left_alpha(i,2) = (left_ang(i+1,2)-2*left_ang(i,2)+left_ang(i-1,2))/(.01^2);
    left_alpha(i,3) = (left_ang(i+1,3)-2*left_ang(i,3)+left_ang(i-1,3))/(.01^2);
end

% Segment Omega and alpha
right_seg_omega = zeros(length(right_globe2ana),3);
left_seg_omega = zeros(length(right_globe2ana),3);
right_seg_alpha = zeros(length(right_seg_omega),3);
left_seg_alpha = zeros(length(left_seg_omega),3);

% Angular Velocity
for i = 1:length(right_globe2ana)
    right_seg_omega(i,1) = cos(right_ang_beta(i))*cos(right_ang_gamma(i))*right_omega(i,1)+sin(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,2) = -cos(right_ang_beta(i))*sin(right_ang_gamma(i))*right_omega(i,1)+cos(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,3) = sin(right_ang_beta(i))*right_omega(i,1)+right_omega(i,3);
    left_seg_omega(i,1) = cos(left_ang_beta(i))*cos(left_ang_gamma(i))*left_omega(i,1)+sin(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,2) = -cos(left_ang_beta(i))*sin(left_ang_gamma(i))*left_omega(i,1)+cos(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,3) = sin(left_ang_beta(i))*left_omega(i,1)+left_omega(i,3);
end

% Angular Acceleration
for i = 2:length(right_ang_alpha)-1
    right_seg_alpha(i,1) = (right_seg_omega(i+1,1)-right_seg_omega(i-1,1))/.02;
    right_seg_alpha(i,2) = (right_seg_omega(i+1,2)-right_seg_omega(i-1,2))/.02;
    right_seg_alpha(i,3) = (right_seg_omega(i+1,3)-right_seg_omega(i-1,3))/.02;
    left_seg_alpha(i,1) = (left_seg_omega(i+1,1)-left_seg_omega(i-1,1))/.02;
    left_seg_alpha(i,2) = (left_seg_omega(i+1,2)-left_seg_omega(i-1,2))/.02;
    left_seg_alpha(i,3) = (left_seg_omega(i+1,3)-left_seg_omega(i-1,3))/.02;
end

end