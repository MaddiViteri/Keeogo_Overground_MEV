% Test Data
% static_markerStruct = static.neutral_static.markers;
% markerStruct = trial.neutral_control.markers;
%% Static Calibration
function [right_globe2ana,left_globe2ana,right_ana2globe,left_ana2globe,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha,RKJC_ana,LKJC_ana] = calibrationmatrix_thigh(static_markerStruct,markerStruct)
%% Marker to anatomical  

% Set Filter
fs_mocap = 100;
    fs_nyq_mo = fs_mocap/2;
    fc_mo_low = 6;
    Wn_mo = fc_mo_low/fs_nyq_mo;
    [cb,ca] = butter(4,Wn_mo,'low');
%Retrive Variables from imported data
    %fix to use with any file
    LASIS = static_markerStruct.LASIS;
    RASIS = static_markerStruct.RASIS;

%determine # of data points for for loop
dataSize = length(LASIS);

%Create matrix of same length as data
Right_HJC = zeros(dataSize,3);
Left_HJC = zeros(dataSize,3);
MidpointBtwnHip = zeros(dataSize,3);

%%
%run for loop
for i = 1:dataSize

    %1) Calculate origin from markers
        origin = (LASIS(i,:) + RASIS(i,:))/2;
        MidpointBtwnHip (i,:) = origin;

    %2) Find ASIS distance from RASIS and LASIS 
        ASIS_Distance = norm(RASIS(i,:)-LASIS(i,:));
        
    %3) RHJC/LHJC calculations
        Right_HJC (i,:) = [-0.36*ASIS_Distance+origin(1),-0.19*ASIS_Distance+origin(2), -0.3*ASIS_Distance+origin(3)];
        Left_HJC (i,:) = [0.36*ASIS_Distance+origin(1),-0.19*ASIS_Distance+origin(2), -0.3*ASIS_Distance+origin(3)];
        
end

hjc.left = mean(Left_HJC,1);
hjc.right = mean(Right_HJC,1);

% Knee Joint Center
kjc.left(:,1) = mean([mean(static_markerStruct.LLatEpi(:,1),"omitnan"),...
    mean(static_markerStruct.LMedEpi(:,1),"omitnan")],2);
kjc.left(:,2) = mean([mean(static_markerStruct.LLatEpi(:,2),"omitnan"),...
    mean(static_markerStruct.LMedEpi(:,2),"omitnan")],2);
kjc.left(:,3) = mean([mean(static_markerStruct.LLatEpi(:,3),"omitnan"),...
    mean(static_markerStruct.LMedEpi(:,3),"omitnan")],2);
kjc.right(:,1) = mean([mean(static_markerStruct.RLatEpi(:,1),"omitnan"),...
    mean(static_markerStruct.RMedEpi(:,1),"omitnan")],2);
kjc.right(:,2) = mean([mean(static_markerStruct.RLatEpi(:,2),"omitnan"),...
    mean(static_markerStruct.RMedEpi(:,2),"omitnan")],2);
kjc.right(:,3) = mean([mean(static_markerStruct.RLatEpi(:,3),"omitnan"),...
    mean(static_markerStruct.RMedEpi(:,3),"omitnan")],2);

% Static COM position
right.static_COM = [0.567*hjc.right(1)+0.433*kjc.right(1),0.567*hjc.right(2)+0.433*kjc.right(2),...
    0.567*hjc.right(3)+0.433*kjc.right(3)];
left.static_COM = [0.567*hjc.left(1)+0.433*kjc.left(1),0.567*hjc.left(2)+0.433*kjc.left(2),...
    0.567*hjc.left(3)+0.433*kjc.left(3)];

%% Segment Coordinate System
% Right 
right_thigh_xan = (mean(static_markerStruct.RMedEpi,"omitnan")-mean(static_markerStruct.RLatEpi,"omitnan"))/...
    norm((mean(static_markerStruct.RMedEpi,"omitnan")-mean(static_markerStruct.RLatEpi,"omitnan")));
right_thigh_yan = [hjc.right(1)-kjc.right(1),hjc.right(2)-kjc.right(2),hjc.right(3)-kjc.right(3)];
right_thigh_zan = cross(right_thigh_xan,right_thigh_yan);
right_thigh_xan = cross(right_thigh_yan,right_thigh_zan);
right_thigh_xan = right_thigh_xan/norm(right_thigh_xan);
right_thigh_yan = right_thigh_yan/norm(right_thigh_yan);
right_thigh_zan = right_thigh_zan/norm(right_thigh_zan);

% Left
left_thigh_xan = (mean(static_markerStruct.LLatEpi,"omitnan")-mean(static_markerStruct.LMedEpi,"omitnan"))/...
    norm((mean(static_markerStruct.LLatEpi,"omitnan")-mean(static_markerStruct.LMedEpi,"omitnan")));
left_thigh_yan = [hjc.left(1)-kjc.left(1),hjc.left(2)-kjc.left(2),hjc.left(3)-kjc.left(3)];
left_thigh_zan = cross(left_thigh_xan,left_thigh_yan);
left_thigh_xan = cross(left_thigh_yan,left_thigh_zan);
left_thigh_xan = left_thigh_xan/norm(left_thigh_xan);
left_thigh_yan = left_thigh_yan/norm(left_thigh_yan);
left_thigh_zan = left_thigh_zan/norm(left_thigh_zan);

% Matrices
RLG2A = [right_thigh_xan;right_thigh_yan;right_thigh_zan];
LLG2A = [left_thigh_xan;left_thigh_yan;left_thigh_zan];
RLA2G = RLG2A';
LLA2G = LLG2A';

%% Marker position relative to anatomical vector position
%% HJC
% Right
right_thigh2hip_xal = hjc.right(1) - right.static_COM(1);
right_thigh2hip_yal = hjc.right(2) - right.static_COM(2);
right_thigh2hip_zal = hjc.right(3) - right.static_COM(3);
right_thigh2hip_al = [right_thigh2hip_xal;right_thigh2hip_yal;right_thigh2hip_zal];
right_an_hip = RLG2A*right_thigh2hip_al;

% Left
left_thigh2hip_xal = hjc.left(1) - left.static_COM(1);
left_thigh2hip_yal = hjc.left(2) - left.static_COM(2);
left_thigh2hip_zal = hjc.left(3) - left.static_COM(3);
left_thigh2hip_al = [left_thigh2hip_xal;left_thigh2hip_yal;left_thigh2hip_zal];
left_an_hip = LLG2A*left_thigh2hip_al;

%% KJC
% Right
right_thigh2kjc_xal = kjc.right(1) - right.static_COM(1);
right_thigh2kjc_yal = kjc.right(2) - right.static_COM(2);
right_thigh2kjc_zal = kjc.right(3) - right.static_COM(3);
right_thigh2kjc_al = [right_thigh2kjc_xal;right_thigh2kjc_yal;right_thigh2kjc_zal];
right_an_kjc = RLG2A*right_thigh2kjc_al;

% Left
left_thigh2kjc_xal = kjc.left(1) - left.static_COM(1);
left_thigh2kjc_yal = kjc.left(2) - left.static_COM(2);
left_thigh2kjc_zal = kjc.left(3) - left.static_COM(3);
left_thigh2kjc_al = [left_thigh2kjc_xal;left_thigh2kjc_yal;left_thigh2kjc_zal];
left_an_kjc = LLG2A*left_thigh2kjc_al;
%% Cluster center
% Right
right_thigh_cluster = mean([static_markerStruct.RAntLatThigh;static_markerStruct.RAntMedThigh;...
    static_markerStruct.RPosLatThigh;static_markerStruct.RPosMedThigh],"omitnan");
right_thigh2cluster_xal = mean(right_thigh_cluster(:,1),"omitnan") - right.static_COM(1);
right_thigh2cluster_yal = mean(right_thigh_cluster(:,2),"omitnan") - right.static_COM(2);
right_thigh2cluster_zal = mean(right_thigh_cluster(:,3),"omitnan") - right.static_COM(3);
right_thigh2cluster_al = [right_thigh2cluster_xal;right_thigh2cluster_yal;right_thigh2cluster_zal];
right_an_cluster = RLG2A*right_thigh2cluster_al;

% Left
left_thigh_cluster = mean([static_markerStruct.LAntLatThigh;static_markerStruct.LAntMedThigh;...
    static_markerStruct.LPosLatThigh;static_markerStruct.LPosMedThigh],"omitnan");
left_thigh2cluster_xal = mean(left_thigh_cluster(:,1),"omitnan") - left.static_COM(1);
left_thigh2cluster_yal = mean(left_thigh_cluster(:,2),"omitnan") - left.static_COM(2);
left_thigh2cluster_zal = mean(left_thigh_cluster(:,3),"omitnan") - left.static_COM(3);
left_thigh2cluster_al = [left_thigh2cluster_xal;left_thigh2cluster_yal;left_thigh2cluster_zal];
left_an_cluster = LLG2A*left_thigh2cluster_al;

%% Ant Lat Thigh
% Right
right_thigh2antlat_xal = mean(static_markerStruct.RAntLatThigh(:,1),"omitnan") - right.static_COM(1);
right_thigh2antlat_yal = mean(static_markerStruct.RAntLatThigh(:,2),"omitnan") - right.static_COM(2);
right_thigh2antlat_zal = mean(static_markerStruct.RAntLatThigh(:,3),"omitnan") - right.static_COM(3);
right_thigh2antlat_al = [right_thigh2antlat_xal;right_thigh2antlat_yal;right_thigh2antlat_zal];
right_an_antlat = RLG2A*right_thigh2antlat_al;

% Left
left_thigh2antlat_xal = mean(static_markerStruct.LAntLatThigh(:,1),"omitnan") - left.static_COM(1);
left_thigh2antlat_yal = mean(static_markerStruct.LAntLatThigh(:,2),"omitnan") - left.static_COM(2);
left_thigh2antlat_zal = mean(static_markerStruct.LAntLatThigh(:,3),"omitnan") - left.static_COM(3);
left_thigh2antlat_al = [left_thigh2antlat_xal;left_thigh2antlat_yal;left_thigh2antlat_zal];
left_an_antlat = LLG2A*left_thigh2antlat_al;

%% Ant med thigh
% Right
right_thigh2antmed_xal = mean(static_markerStruct.RAntMedThigh(:,1),"omitnan") - right.static_COM(1);
right_thigh2antmed_yal = mean(static_markerStruct.RAntMedThigh(:,2),"omitnan") - right.static_COM(2);
right_thigh2antmed_zal = mean(static_markerStruct.RAntMedThigh(:,3),"omitnan") - right.static_COM(3);
right_thigh2antmed_al = [right_thigh2antmed_xal;right_thigh2antmed_yal;right_thigh2antmed_zal];
right_an_antmed = RLG2A*right_thigh2antmed_al;

% Left
left_thigh2antmed_xal = mean(static_markerStruct.LAntMedThigh(:,1),"omitnan") - left.static_COM(1);
left_thigh2antmed_yal = mean(static_markerStruct.LAntMedThigh(:,2),"omitnan") - left.static_COM(2);
left_thigh2antmed_zal = mean(static_markerStruct.LAntMedThigh(:,3),"omitnan") - left.static_COM(3);
left_thigh2antmed_al = [left_thigh2antmed_xal;left_thigh2antmed_yal;left_thigh2antmed_zal];
left_an_antmed = LLG2A*left_thigh2antmed_al;

%% Pos Med thigh
% Right
right_thigh2posmed_xal = mean(static_markerStruct.RPosMedThigh(:,1),"omitnan") - right.static_COM(1);
right_thigh2posmed_yal = mean(static_markerStruct.RPosMedThigh(:,2),"omitnan") - right.static_COM(2);
right_thigh2posmed_zal = mean(static_markerStruct.RPosMedThigh(:,3),"omitnan") - right.static_COM(3);
right_thigh2posmed_al = [right_thigh2posmed_xal;right_thigh2posmed_yal;right_thigh2posmed_zal];
right_an_posmed = RLG2A*right_thigh2posmed_al;

% Left
left_thigh2posmed_xal = mean(static_markerStruct.LPosMedThigh(:,1),"omitnan") - left.static_COM(1);
left_thigh2posmed_yal = mean(static_markerStruct.LPosMedThigh(:,2),"omitnan") - left.static_COM(2);
left_thigh2posmed_zal = mean(static_markerStruct.LPosMedThigh(:,3),"omitnan") - left.static_COM(3);
left_thigh2posmed_al = [left_thigh2posmed_xal;left_thigh2posmed_yal;left_thigh2posmed_zal];
left_an_posmed = LLG2A*left_thigh2posmed_al;

%% Marker Matrix Coordinate System
% Right
right_local_zm = (right_an_antmed - right_an_antlat)/norm(right_an_antmed - right_an_antlat);
right_local_a = (right_an_posmed - right_an_antlat)/norm(right_an_posmed - right_an_antlat);
right_local_ym = cross(right_local_zm,right_local_a);
right_local_xm = cross(right_local_ym,right_local_zm);
right_local_xm = right_local_xm/norm(right_local_xm);
right_local_ym = right_local_ym/norm(right_local_ym);
right_local_zm = right_local_zm/norm(right_local_zm);
right_LM2A = [right_local_xm,right_local_ym,right_local_zm];

% Left
left_local_zm = (left_an_antmed - left_an_antlat)/norm(left_an_antmed - left_an_antlat);
left_local_a = (left_an_posmed - left_an_antlat)/norm(left_an_posmed - left_an_antlat);
left_local_ym = cross(left_local_zm,left_local_a);
left_local_xm = cross(left_local_ym,left_local_zm);
left_local_xm = left_local_xm/norm(left_local_xm);
left_local_ym = left_local_ym/norm(left_local_ym);
left_local_zm = left_local_zm/norm(left_local_zm);
left_LM2A = [left_local_xm,left_local_ym,left_local_zm];

%% Global to Marker
% Create matrix to get from the global to the marker coordinate system
% during the actual trial

% Initialize Data
right_globe_am = zeros(length(markerStruct.RLatMal),3);
right_globe_xm = zeros(length(markerStruct.RLatMal),3);
right_globe_ym = zeros(length(markerStruct.RLatMal),3);
right_globe_zm = zeros(length(markerStruct.RLatMal),3);
left_globe_am = zeros(length(markerStruct.RLatMal),3);
left_globe_xm = zeros(length(markerStruct.RLatMal),3);
left_globe_ym = zeros(length(markerStruct.RLatMal),3);
left_globe_zm = zeros(length(markerStruct.RLatMal),3);
right_globe2mark = zeros(3,3,length(markerStruct.RLatMal));
left_globe2mark = zeros(3,3,length(markerStruct.RLatMal));

right_globe_zm = (markerStruct.RAntMedThigh - markerStruct.RAntLatThigh)./...
    sqrt(sum((markerStruct.RAntMedThigh - markerStruct.RAntLatThigh).^2,2));
right_globe_am = (markerStruct.RPosMedThigh - markerStruct.RAntLatThigh)./...
    sqrt(sum((markerStruct.RPosMedThigh - markerStruct.RAntLatThigh).^2,2));
right_globe_ym = cross(right_globe_zm,right_globe_am,2);
right_globe_ym = right_globe_ym./sqrt(sum(right_globe_ym.^2,2));
right_globe_xm = cross(right_globe_ym, right_globe_zm,2);
right_globe_xm = right_globe_xm./sqrt(sum(right_globe_xm.^2,2));

left_globe_zm = (markerStruct.LAntMedThigh - markerStruct.LAntLatThigh)./...
    sqrt(sum((markerStruct.LAntMedThigh - markerStruct.LAntLatThigh).^2,2));
left_globe_am = (markerStruct.LPosMedThigh - markerStruct.LAntLatThigh)./...
    sqrt(sum((markerStruct.LPosMedThigh - markerStruct.LAntLatThigh).^2,2));
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

%% Anatomical to Global
for i = 1:length(right_globe2mark)
    right_ana2globe(:,:,i) = right_globe2ana(:,:,i)';
    left_ana2globe(:,:,i) = left_globe2ana(:,:,i)';
end
%% Angle Calculation
% General premise of calculating the joint/segment angles from rotation
% matrix
% syms x y z
% mat1 = [1 0 0;0 cos(x) sin(x); 0 -sin(x) cos(x)];
% mat2 = [cos(y) 0 -sin(y);0 1 0;-sin(y) 0 cos(y)];
% mat3 = [cos(z) sin(z) 0;-sin(z) cos(z) 0;0 0 1];
% rot = mat3*mat2*mat1;

%% Calculate joint angles
right_ang_alpha = unwrap(atan2(-right_globe2ana(3,2,:),right_globe2ana(3,3,:)));
right_ang_beta = unwrap(atan2(right_globe2ana(3,1,:),sqrt((right_globe2ana(1,1,:)).^2+(right_globe2ana(2,1,:)).^2)));
right_ang_gamma = unwrap(atan2(-right_globe2ana(2,1,:),right_globe2ana(1,1,:)));
right_ang = [right_ang_alpha(:),right_ang_beta(:),right_ang_gamma(:)];
left_ang_alpha = unwrap(atan2(-left_globe2ana(3,2,:),left_globe2ana(3,3,:)));
left_ang_beta = unwrap(atan2(left_globe2ana(3,1,:),sqrt((left_globe2ana(1,1,:)).^2+(left_globe2ana(2,1,:)).^2)));
left_ang_gamma = unwrap(atan2(-left_globe2ana(2,1,:),left_globe2ana(1,1,:)));
left_ang = [left_ang_alpha(:),left_ang_beta(:),left_ang_gamma(:)];

%% Filter Data
right_ang(isnan(right_ang)) = 0;
right_ang = filtfilt(cb,ca,right_ang);
left_ang(isnan(left_ang)) = 0;
left_ang = filtfilt(cb,ca,left_ang);

%% Initialize variables
right_omega = zeros(length(right_ang_alpha),3);
left_omega = zeros(length(left_ang_alpha),3);
right_alpha = zeros(length(right_ang_alpha),3);
left_alpha = zeros(length(left_ang_alpha),3);
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
for i = 1:length(right_globe2ana)
    right_seg_omega(i,1) = cos(right_ang_beta(i))*cos(right_ang_gamma(i))*right_omega(i,1)+sin(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,2) = -cos(right_ang_beta(i))*sin(right_ang_gamma(i))*right_omega(i,1)+cos(right_ang_gamma(i))*right_omega(i,2);
    right_seg_omega(i,3) = sin(right_ang_beta(i))*right_omega(i,1)+right_omega(i,3);
    left_seg_omega(i,1) = cos(left_ang_beta(i))*cos(left_ang_gamma(i))*left_omega(i,1)+sin(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,2) = -cos(left_ang_beta(i))*sin(left_ang_gamma(i))*left_omega(i,1)+cos(left_ang_gamma(i))*left_omega(i,2);
    left_seg_omega(i,3) = sin(left_ang_beta(i))*left_omega(i,1)+left_omega(i,3);
end
right_seg_alpha = zeros(length(right_seg_omega),3);
left_seg_alpha = zeros(length(left_seg_omega),3);

% Segment angular acceleration in segment coordinate system
for i = 2:length(right_ang_alpha)-1
    right_seg_alpha(i,1) = (right_seg_omega(i+1,1)-right_seg_omega(i-1,1))/.02;
    right_seg_alpha(i,2) = (right_seg_omega(i+1,2)-right_seg_omega(i-1,2))/.02;
    right_seg_alpha(i,3) = (right_seg_omega(i+1,3)-right_seg_omega(i-1,3))/.02;
    left_seg_alpha(i,1) = (left_seg_omega(i+1,1)-left_seg_omega(i-1,1))/.02;
    left_seg_alpha(i,2) = (left_seg_omega(i+1,2)-left_seg_omega(i-1,2))/.02;
    left_seg_alpha(i,3) = (left_seg_omega(i+1,3)-left_seg_omega(i-1,3))/.02;
end

%% Cluster position
% Right
right_thigh_cluster_globe(:,1) = mean([markerStruct.RAntLatThigh(:,1),markerStruct.RAntMedThigh(:,1), ...
    markerStruct.RPosLatThigh(:,1),markerStruct.RPosMedThigh(:,1)],2,"omitnan");
right_thigh_cluster_globe(:,2) = mean([markerStruct.RAntLatThigh(:,2),markerStruct.RAntMedThigh(:,2), ...
    markerStruct.RPosLatThigh(:,2),markerStruct.RPosMedThigh(:,2)],2,"omitnan");
right_thigh_cluster_globe(:,3) = mean([markerStruct.RAntLatThigh(:,3),markerStruct.RAntMedThigh(:,3), ...
    markerStruct.RPosLatThigh(:,3),markerStruct.RPosMedThigh(:,3)],2,"omitnan");

% Left
left_thigh_cluster_globe(:,1) = mean([markerStruct.LAntLatThigh(:,1),markerStruct.LAntMedThigh(:,1), ...
    markerStruct.LPosLatThigh(:,1),markerStruct.LPosMedThigh(:,1)],2,"omitnan");
left_thigh_cluster_globe(:,2) = mean([markerStruct.LAntLatThigh(:,2),markerStruct.LAntMedThigh(:,2), ...
    markerStruct.LPosLatThigh(:,2),markerStruct.LPosMedThigh(:,2)],2,"omitnan");
left_thigh_cluster_globe(:,3) = mean([markerStruct.LAntLatThigh(:,3),markerStruct.LAntMedThigh(:,3), ...
    markerStruct.LPosLatThigh(:,3),markerStruct.LPosMedThigh(:,3)],2,"omitnan");

%% Knee joint position
for i = 1:length(right_ana2globe)
    RKJC_ana(i,:) = ((right_ana2globe(:,:,i)*(-right_an_cluster))' +...
        (right_ana2globe(:,:,i)*right_an_kjc)' +right_thigh_cluster_globe(i,:))/1000;
    LKJC_ana(i,:) = ((left_ana2globe(:,:,i)*(-left_an_cluster))' +...
        (left_ana2globe(:,:,i)*left_an_kjc)' +left_thigh_cluster_globe(i,:))/1000;
end

end