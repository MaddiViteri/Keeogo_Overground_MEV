% Test Data
% static_markerStruct = static.neutral_static.markers;
% markerStruct = trial.neutral_control.markers;
%% Static Calibration
% Load OpenSim libs
function [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha] = calibrationmatrix_6_27(static_markerStruct,markerStruct)
%% Marker to anatomical

% Set filters
fs_mocap = 100;
    fs_nyq_mo = fs_mocap/2;
    fc_mo_low = 6;
    Wn_mo = fc_mo_low/fs_nyq_mo;
    [cb,ca] = butter(4,Wn_mo,'low');

% Determine Joint Centers
ajc.left(:,1) = mean([mean(static_markerStruct.LLatMal(:,1),"omitnan"),mean(static_markerStruct.LMedMal(:,1),"omitnan")],2);
ajc.left(:,2) = mean([mean(static_markerStruct.LLatMal(:,2),"omitnan"),mean(static_markerStruct.LMedMal(:,2),"omitnan")]);
ajc.left(:,3) = mean([mean(static_markerStruct.LLatMal(:,3),"omitnan"),mean(static_markerStruct.LMedMal(:,3),"omitnan")]);
ajc.right(:,1) = mean([mean(static_markerStruct.RLatMal(:,1),"omitnan"),mean(static_markerStruct.RMedMal(:,1),"omitnan")]);
ajc.right(:,2) = mean([mean(static_markerStruct.RLatMal(:,2),"omitnan"),mean(static_markerStruct.RMedMal(:,2),"omitnan")]);
ajc.right(:,3) = mean([mean(static_markerStruct.RLatMal(:,3),"omitnan"),mean(static_markerStruct.RMedMal(:,3),"omitnan")]);
kjc.left(:,1) = mean([mean(static_markerStruct.LLatEpi(:,1),"omitnan"),mean(static_markerStruct.LMedEpi(:,1),"omitnan")]);
kjc.left(:,2) = mean([mean(static_markerStruct.LLatEpi(:,2),"omitnan"),mean(static_markerStruct.LMedEpi(:,2),"omitnan")]);
kjc.left(:,3) = mean([mean(static_markerStruct.LLatEpi(:,3),"omitnan"),mean(static_markerStruct.LMedEpi(:,3),"omitnan")]);
kjc.right(:,1) = mean([mean(static_markerStruct.RLatEpi(:,1),"omitnan"),mean(static_markerStruct.RMedEpi(:,1),"omitnan")]);
kjc.right(:,2) = mean([mean(static_markerStruct.RLatEpi(:,2),"omitnan"),mean(static_markerStruct.RMedEpi(:,2),"omitnan")]);
kjc.right(:,3) = mean([mean(static_markerStruct.RLatEpi(:,3),"omitnan"),mean(static_markerStruct.RMedEpi(:,3),"omitnan")]);

% Determine center of mass from static trial 
right.static_COM = [0.567*kjc.right(1)+0.433*ajc.right(1),0.567*kjc.right(2)+0.433*ajc.right(2),...
    0.567*kjc.right(3)+0.433*ajc.right(3)];
left.static_COM = [0.567*kjc.left(1)+0.433*ajc.left(1),0.567*kjc.left(2)+0.433*ajc.left(2),...
    0.567*kjc.left(3)+0.433*ajc.left(3)];
% Anatomical (an) Axes calculation
% Right
right_shank_xan = (mean(static_markerStruct.RMedMal,"omitnan")-mean(static_markerStruct.RLatMal,"omitnan"))/...
    norm((mean(static_markerStruct.RMedMal,"omitnan")-mean(static_markerStruct.RLatMal,"omitnan")));
right_shank_yan = [kjc.right(1)-ajc.right(1),kjc.right(2)-ajc.right(2),kjc.right(3)-ajc.right(3)];
right_shank_zan = cross(right_shank_xan,right_shank_yan);
right_shank_xan = cross(right_shank_yan,right_shank_zan);
right_shank_xan = right_shank_xan/norm(right_shank_xan);
right_shank_yan = right_shank_yan/norm(right_shank_yan);
right_shank_zan = right_shank_zan/norm(right_shank_zan);

% Left
left_shank_xan = (mean(static_markerStruct.LLatMal,"omitnan")-mean(static_markerStruct.LMedMal,"omitnan"))/...
    norm((mean(static_markerStruct.LLatMal,"omitnan")-mean(static_markerStruct.LMedMal,"omitnan")));
left_shank_yan = [kjc.left(1)-ajc.left(1),kjc.left(2)-ajc.left(2),kjc.left(3)-ajc.left(3)];
left_shank_zan = cross(left_shank_xan,left_shank_yan);
left_shank_xan = cross(left_shank_yan,left_shank_zan);
left_shank_xan = left_shank_xan/norm(left_shank_xan);
left_shank_yan = left_shank_yan/norm(left_shank_yan);
left_shank_zan = left_shank_zan/norm(left_shank_zan);

% Matrices
RLG2A = [right_shank_xan;right_shank_yan;right_shank_zan];
LLG2A = [left_shank_xan;left_shank_yan;left_shank_zan];
RLA2G = RLG2A';
LLA2G = LLG2A';
% Marker's relative anatomical vector position
% Ankle
right_shank2ank_xal = ajc.right(1) - right.static_COM(1);
right_shank2ank_yal = ajc.right(2) - right.static_COM(2);
right_shank2ank_zal = ajc.right(3) - right.static_COM(3);
right_shank2ank_al = [right_shank2ank_xal;right_shank2ank_yal;right_shank2ank_zal];
right_an_ank = RLG2A*right_shank2ank_al;
left_shank2ank_xal = ajc.left(1) - left.static_COM(1);
left_shank2ank_yal = ajc.left(2) - left.static_COM(2);
left_shank2ank_zal = ajc.left(3) - left.static_COM(3);
left_shank2ank_al = [left_shank2ank_xal;left_shank2ank_yal;left_shank2ank_zal];
left_an_ank = LLG2A*left_shank2ank_al;

% Lateral Malleolus
right_shank2latmal_xal = mean(static_markerStruct.RLatMal(:,1),"omitnan") - right.static_COM(1);
right_shank2latmal_yal = mean(static_markerStruct.RLatMal(:,2),"omitnan") - right.static_COM(2);
right_shank2latmal_zal = mean(static_markerStruct.RLatMal(:,3),"omitnan") - right.static_COM(3);
right_shank2latmal_al = [right_shank2latmal_xal;right_shank2latmal_yal;right_shank2latmal_zal];
right_an_latmal = RLG2A*right_shank2latmal_al;
left_shank2latmal_xal = mean(static_markerStruct.LLatMal(:,1),"omitnan") - left.static_COM(1);
left_shank2latmal_yal = mean(static_markerStruct.LLatMal(:,2),"omitnan") - left.static_COM(2);
left_shank2latmal_zal = mean(static_markerStruct.LLatMal(:,3),"omitnan") - left.static_COM(3);
left_shank2latmal_al = [left_shank2latmal_xal;left_shank2latmal_yal;left_shank2latmal_zal];
left_an_latmal = LLG2A*left_shank2latmal_al;

% Medial Upper Shank
right_shank2medepi_xal = mean(static_markerStruct.RSupMedShank(:,1),"omitnan") - right.static_COM(1);
right_shank2medepi_yal = mean(static_markerStruct.RSupMedShank(:,2),"omitnan") - right.static_COM(2);
right_shank2medepi_zal = mean(static_markerStruct.RSupMedShank(:,3),"omitnan") - right.static_COM(3);
right_shank2medepi_al = [right_shank2medepi_xal;right_shank2medepi_yal;right_shank2medepi_zal];
right_an_medepi = RLG2A*right_shank2medepi_al;
left_shank2medepi_xal = mean(static_markerStruct.LSupMedShank(:,1),"omitnan") - left.static_COM(1);
left_shank2medepi_yal = mean(static_markerStruct.LSupMedShank(:,2),"omitnan") - left.static_COM(2);
left_shank2medepi_zal = mean(static_markerStruct.LSupMedShank(:,3),"omitnan") - left.static_COM(3);
left_shank2medepi_al = [left_shank2medepi_xal;left_shank2medepi_yal;left_shank2medepi_zal];
left_an_medepi = LLG2A*left_shank2medepi_al;

% Medial Lower Shank
right_shank2LLShank_xal = mean(static_markerStruct.RInfMedShank(:,1),"omitnan") - right.static_COM(1);
right_shank2LLShank_yal = mean(static_markerStruct.RInfMedShank(:,2),"omitnan") - right.static_COM(2);
right_shank2LLShank_zal = mean(static_markerStruct.RInfMedShank(:,3),"omitnan") - right.static_COM(3);
right_shank2LLShank_al = [right_shank2LLShank_xal;right_shank2LLShank_yal;right_shank2LLShank_zal];
right_an_LLShank = RLG2A*right_shank2LLShank_al;
left_shank2LRShank_xal = mean(static_markerStruct.LInfMedShank(:,1),"omitnan") - left.static_COM(1);
left_shank2LRShank_yal = mean(static_markerStruct.LInfMedShank(:,2),"omitnan") - left.static_COM(2);
left_shank2LRShank_zal = mean(static_markerStruct.LInfMedShank(:,3),"omitnan") - left.static_COM(3);
left_shank2LRShank_al = [left_shank2LRShank_xal;left_shank2LRShank_yal;left_shank2LRShank_zal];
left_an_LRShank = LLG2A*left_shank2LRShank_al;

% Marker Matrix
% Right
right_local_zm = (right_an_medepi - right_an_latmal)/norm(right_an_medepi - right_an_latmal);
right_local_a = (right_an_LLShank - right_an_latmal)/norm(right_an_LLShank - right_an_latmal);
right_local_ym = cross(right_local_zm,right_local_a);
right_local_xm = cross(right_local_ym,right_local_zm);
right_local_xm = right_local_xm/norm(right_local_xm);
right_local_ym = right_local_ym/norm(right_local_ym);
right_local_zm = right_local_zm/norm(right_local_zm);
right_LM2A = [right_local_xm,right_local_ym,right_local_zm];
% Left
left_local_zm = (left_an_medepi - left_an_latmal)/norm(left_an_medepi - left_an_latmal);
left_local_a = (left_an_LRShank - left_an_latmal)/norm(left_an_LRShank - left_an_latmal);
left_local_ym = cross(left_local_zm,left_local_a);
left_local_xm = cross(left_local_ym,left_local_zm);
left_local_xm = left_local_xm/norm(left_local_xm);
left_local_ym = left_local_ym/norm(left_local_ym);
left_local_zm = left_local_zm/norm(left_local_zm);
left_LM2A = [left_local_xm,left_local_ym,left_local_zm];
%% Global to Marker
% Create rotation matrix to get from global positioning to marker position
% during the actual trial as mediary step

% Initialize Matrices
right_globe2mark = zeros(3,3,length(markerStruct.RLatMal));
left_globe2mark = zeros(3,3,length(markerStruct.RLatMal));

% Matrices to Global to Marker for Trial
right_globe_zm = (markerStruct.RSupMedShank - markerStruct.RLatMal)./...
    sqrt(sum((markerStruct.RSupMedShank - markerStruct.RLatMal).^2,2));
right_globe_am = (markerStruct.RInfMedShank - markerStruct.RLatMal)./...
    sqrt(sum((markerStruct.RInfMedShank - markerStruct.RLatMal).^2,2));
right_globe_ym = cross(right_globe_zm,right_globe_am,2);
right_globe_ym = right_globe_ym./sqrt(sum(right_globe_ym.^2,2));
right_globe_xm = cross(right_globe_ym, right_globe_zm,2);
right_globe_xm = right_globe_xm./sqrt(sum(right_globe_xm.^2,2));

left_globe_zm = (markerStruct.LSupMedShank - markerStruct.LLatMal)./...
    sqrt(sum((markerStruct.LSupMedShank - markerStruct.LLatMal).^2,2));
left_globe_am = (markerStruct.LInfMedShank - markerStruct.LLatMal)./...
    sqrt(sum((markerStruct.LInfMedShank - markerStruct.LLatMal).^2,2));
left_globe_ym = cross(left_globe_zm,left_globe_am,2);
left_globe_ym = left_globe_ym./sqrt(sum(left_globe_ym.^2,2));
left_globe_xm = cross(left_globe_ym,left_globe_zm,2);
left_globe_xm = left_globe_xm./sqrt(sum(left_globe_xm.^2,2));

for i = 1:length(left_globe_zm)
    right_globe2mark(:,:,i) =  [right_globe_xm(i,:);right_globe_ym(i,:);right_globe_zm(i,:)];
    left_globe2mark(:,:,i) =  [left_globe_xm(i,:);left_globe_ym(i,:);left_globe_zm(i,:)];
end
%% Global to Anatomical
for i = 1:length(right_globe2mark)
    right_globe2ana(:,:,i) = right_LM2A*right_globe2mark(:,:,i);
    left_globe2ana(:,:,i) = left_LM2A*left_globe2mark(:,:,i);
end

%% Angle Calculation
% syms x y z
% mat1 = [1 0 0;0 cos(x) sin(x); 0 -sin(x) cos(x)];
% mat2 = [cos(y) 0 -sin(y);0 1 0;-sin(y) 0 cos(y)];
% mat3 = [cos(z) sin(z) 0;-sin(z) cos(z) 0;0 0 1];
% rot = mat3*mat2*mat1;

right_ang_alpha = unwrap(atan2(-right_globe2ana(3,2,:),right_globe2ana(3,3,:)));
right_ang_beta = unwrap(atan2(right_globe2ana(3,1,:),sqrt((right_globe2ana(1,1,:)).^2+(right_globe2ana(2,1,:)).^2)));
right_ang_gamma = unwrap(atan2(-right_globe2ana(2,1,:),right_globe2ana(1,1,:)));
right_ang = [right_ang_alpha(:),right_ang_beta(:),right_ang_gamma(:)];
left_ang_alpha = unwrap(atan2(-left_globe2ana(3,2,:),left_globe2ana(3,3,:)));
left_ang_beta = unwrap(atan2(left_globe2ana(3,1,:),sqrt((left_globe2ana(1,1,:)).^2+(left_globe2ana(2,1,:)).^2)));
left_ang_gamma = unwrap(atan2(-left_globe2ana(2,1,:),left_globe2ana(1,1,:)));
left_ang = [left_ang_alpha(:),left_ang_beta(:),left_ang_gamma(:)];

% Filter angle data
right_ang(isnan(right_ang)) = 0;
right_ang = filtfilt(cb,ca,right_ang);
left_ang(isnan(left_ang)) = 0;
left_ang = filtfilt(cb,ca,left_ang);

% Initialize Variables
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

%% Calculate Segment Omega and Alpha
% Segment Angular Velocity
right_seg_omega = zeros(length(right_globe2ana),3);
left_seg_omega = zeros(length(right_globe2ana),3);
for i = 1:length(right_globe2ana)
    right_seg_omega(i,1) = cos(right_ang_beta(i))*cos(right_ang_gamma(i))*right_omega(i,1)+sin(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,2) = -cos(right_ang_beta(i))*sin(right_ang_gamma(i))*right_omega(i,1)+cos(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,3) = sin(right_ang_beta(i))*right_omega(i,1)+right_omega(i,3);
    left_seg_omega(i,1) = cos(left_ang_beta(i))*cos(left_ang_gamma(i))*left_omega(i,1)+sin(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,2) = -cos(left_ang_beta(i))*sin(left_ang_gamma(i))*left_omega(i,1)+cos(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,3) = sin(left_ang_beta(i))*left_omega(i,1)+left_omega(i,3);
end

% Segment Angular Acceleration
right_seg_alpha = zeros(length(right_seg_omega),3);
left_seg_alpha = zeros(length(left_seg_omega),3);
for i = 2:length(right_ang_alpha)-1
    right_seg_alpha(i,1) = (right_seg_omega(i+1,1)-right_seg_omega(i-1,1))/.02;
    right_seg_alpha(i,2) = (right_seg_omega(i+1,2)-right_seg_omega(i-1,2))/.02;
    right_seg_alpha(i,3) = (right_seg_omega(i+1,3)-right_seg_omega(i-1,3))/.02;
    left_seg_alpha(i,1) = (left_seg_omega(i+1,1)-left_seg_omega(i-1,1))/.02;
    left_seg_alpha(i,2) = (left_seg_omega(i+1,2)-left_seg_omega(i-1,2))/.02;
    left_seg_alpha(i,3) = (left_seg_omega(i+1,3)-left_seg_omega(i-1,3))/.02;
end

end