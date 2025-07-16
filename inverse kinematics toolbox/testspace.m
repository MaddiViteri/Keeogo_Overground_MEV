%inverse kinematics toolbox
%Author: Maddi Viteri
%Last Updated: 06/10/2025 (MATLAB R2025a)
%test space for functions

import org.opensim.modeling.*


% pm = PatientManager();
% 
% p1 = PatientData(1, 215, 66, 63, 'male', 'right');
% pm.addPatient(p1);


[filename, filepath] = uigetfile('*.*'); % '*.*' allows all files types to be selected

%[platform] = CARENtxt_to_Struct (filepath, filename); 

[STATICMARKERS, FORCES] = c3d_to_Struct (filepath, filename);

[filename, filepath] = uigetfile('*.*'); % '*.*' allows all files types to be selected

[MARKERS, FORCES] = c3d_to_Struct (filepath, filename);

[FOOT, SHANK, THIGH] = combined_Calibration (STATICMARKERS, MARKERS);