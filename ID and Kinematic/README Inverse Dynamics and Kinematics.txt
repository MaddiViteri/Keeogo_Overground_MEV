
Inverse_3D_Dynamics_7_18.m
- Function calculates the ankle and knee flexion, adduction, rotation moments and the joint forces
- Transforms the kinetic values relative to the the shank coordinate system,

"calibrationmatrix_6_27.m" 
- Uses the Winter 4th Edition to calculate segment coordinate system
- First finds the segments center of mass
- Then determine the segment coordinate system with the calibration markers to find the local global to anatomical matrix
- Next finds the marker position relative to the center of mass
- Use the relative positions to get a local marker to anatomical matrix
- After find the global to marker matrix from the actual trials
- Finally get a get the global to anatomical matrix by multiplying the local marker matrix to anatomical and the global to marker coordinate system
- Code also finds angular position, velocity, and acceleration that is used for inverse dynamics 

"calibrationmatrix_thigh.m"
- Uses Coda Pelvis to determine the hip joint centers
- Also determines knee joint center based on the other things that exist 
- Uses the Winter 4th Edition to calculate segment coordinate system
- First finds the segments center of mass
- Then determine the segment coordinate system with the calibration markers to find the local global to anatomical matrix
- Next finds the marker position relative to the center of mass
- Use the relative positions to get a local marker to anatomical matrix
- After find the global to marker matrix from the actual trials
- Finally get a get the global to anatomical matrix by multiplying the local marker matrix to anatomical and the global to marker coordinate system
- Code also finds angular position, velocity, and acceleration that is used for inverse dynamics 


"calibrationmatrix_foot_6_27.m" 
- Uses the Winter 4th Edition to calculate segment coordinate system
- First finds the segments center of mass
- Then determine the segment coordinate system with the calibration markers to find the local global to anatomical matrix
- Next finds the marker position relative to the center of mass
- Use the relative positions to get a local marker to anatomical matrix
- After find the global to marker matrix from the actual trials
- Finally get a get the global to anatomical matrix by multiplying the local marker matrix to anatomical and the global to marker coordinate system
- Code also finds angular position, velocity, and acceleration that is used for inverse dynamics 

"Step_Cycle_MI_combined_brace_v_keeogo_24_03_13.m"
- Breaks down and normalizes the kinetic measures to 100% of stance phase 

"Step_Length_7_25_23.m"
-Calculates all spatio-temporal measures except for range of motion

"stride_peak_values.m"
- Extracts the peak values of a the measure based on the stage of stance (loading, early stance, mid stnace, and late stance)

"gait_parameters_combined.m"
- Combines the spatiotemporal values of the trials together

"main_parameters_combined.m"
- Combines all the trials together to aggregate the steps and the various values
- Calculates the mean and standarad deviation of them all

"plot_function.m"

- Separate for right and left knee but also combined right and left knee plots
- Plot the mean and standard deviation of the values for the stance phase
