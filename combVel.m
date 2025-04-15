function [data_velocity] = combVel (data_Noraxon)
    
%nss
    tempa = data_Noraxon(1).InsoleL.stepvelocity;
    tempb = data_Noraxon(1).InsoleR.stepvelocity;
    tempc = data_Noraxon(2).InsoleL.stepvelocity;
    tempd = data_Noraxon(2).InsoleR.stepvelocity;

    % Step 1: transpose
    tempa = tempa'; % Transpose to convert 1x16 to 16x1
    tempb = tempb'; % Transpose to convert 1x16 to 16x1
    tempc = tempc'; % Transpose to convert 1x16 to 16x1
    tempd = tempd'; % Transpose to convert 1x16 to 16x1
    
    % Concatenate vertically
    combined_array = [tempa; tempb; tempc; tempd]; 

    % Remove NaN values
    cleaned_array = rmmissing(combined_array);

    % Get unique values
    unique_array = unique(cleaned_array);

    combined_array = unique_array;

    % Step 4: Compute the calculations, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    step_velocities = [data_Noraxon(1).InsoleL.stepvelocity, ...
                  data_Noraxon(1).InsoleR.stepvelocity, ...
                  data_Noraxon(2).InsoleL.stepvelocity, ...
                  data_Noraxon(2).InsoleR.stepvelocity];

    step_endTime = [data_Noraxon(1).InsoleL.stepEndTime, ...
                  data_Noraxon(1).InsoleR.stepEndTime, ...
                  data_Noraxon(2).InsoleL.stepEndTime, ...
                  data_Noraxon(2).InsoleR.stepEndTime];

    %account for outliers
    q1 = median_value - iqr_value;
    q3 = median_value + iqr_value;

    iqr_ = q3 - q1;

    % Step 3: Define outliers as those below (Q1 - 1.5*IQR) or above (Q3 + 1.5*IQR)
    lower_bound = q1 - 1.5 * iqr_;
    upper_bound = q3 + 1.5 * iqr_;
    
    data  = step_velocities;
    % Step 4: Replace outliers with NaN
    data(data < lower_bound | data > upper_bound) = NaN;

    stepVelTime = [data; step_endTime];

    % Find the maximum value and its index
    [max_value, max_index] = max(stepVelTime(1,:));

    median_value = median_value/1000; iqr_value = iqr_value/1000; max_value = max_value /1000;

    %time % shit
    %value from trial 1
    tempStepVel = data_Noraxon(1).InsoleL.stepvelocity, ...
                  data_Noraxon(1).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(1).InsoleL.stepEndTime, ...
                  data_Noraxon(1).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(1).InsoleL.time(data_Noraxon(1).InsoleL.startIndx);
    stop = data_Noraxon(1).InsoleL.time(data_Noraxon(1).InsoleL.endIndx);
    length = stop - start;

    max_tPercent1 = ((max_t-start)/(length))*100;

    %value from trial 2
    tempStepVel = data_Noraxon(2).InsoleL.stepvelocity, ...
                  data_Noraxon(2).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(2).InsoleL.stepEndTime, ...
                  data_Noraxon(2).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(2).InsoleL.time(data_Noraxon(2).InsoleL.startIndx);
    stop = data_Noraxon(2).InsoleL.time(data_Noraxon(2).InsoleL.endIndx);
    length = stop - start;
    
    max_tPercent2 = ((max_t-start)/(length))*100;
            
    %median
    timeMaxVel = median([max_tPercent1,max_tPercent2]);
    
    nss = struct('condition','nss','median',median_value, 'iqr', iqr_value, 'max', max_value, 'timeatpeakvel',timeMaxVel);

%nfs
    tempa = data_Noraxon(3).InsoleL.stepvelocity;
    tempb = data_Noraxon(3).InsoleR.stepvelocity;
    tempc = data_Noraxon(4).InsoleL.stepvelocity;
    tempd = data_Noraxon(4).InsoleR.stepvelocity;

    % Step 1: transpose
    tempa = tempa'; % Transpose to convert 1x16 to 16x1
    tempb = tempb'; % Transpose to convert 1x16 to 16x1
    tempc = tempc'; % Transpose to convert 1x16 to 16x1
    tempd = tempd'; % Transpose to convert 1x16 to 16x1
    
    % Concatenate vertically
    combined_array = [tempa; tempb; tempc; tempd]; 

    % Remove NaN values
    cleaned_array = rmmissing(combined_array);

    % Get unique values
    unique_array = unique(cleaned_array);

    combined_array = unique_array;

    % Step 4: Compute the calculations, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    step_velocities = [data_Noraxon(3).InsoleL.stepvelocity, ...
                  data_Noraxon(3).InsoleR.stepvelocity, ...
                  data_Noraxon(4).InsoleL.stepvelocity, ...
                  data_Noraxon(4).InsoleR.stepvelocity];

    step_endTime = [data_Noraxon(3).InsoleL.stepEndTime, ...
                  data_Noraxon(3).InsoleR.stepEndTime, ...
                  data_Noraxon(4).InsoleL.stepEndTime, ...
                  data_Noraxon(4).InsoleR.stepEndTime];

    %account for outliers
    q1 = median_value - iqr_value;
    q3 = median_value + iqr_value;

    iqr_ = q3 - q1;

    % Step 3: Define outliers as those below (Q1 - 1.5*IQR) or above (Q3 + 1.5*IQR)
    lower_bound = q1 - 1.5 * iqr_;
    upper_bound = q3 + 1.5 * iqr_;
    
    data  = step_velocities;
    % Step 4: Replace outliers with NaN
    data(data < lower_bound | data > upper_bound) = NaN;

    stepVelTime = [data; step_endTime];

    % Find the maximum value and its index
    [max_value, max_index] = max(stepVelTime(1,:));
    
    timeMaxVel = stepVelTime(2,max_index);

    median_value = median_value/1000; iqr_value = iqr_value/1000; max_value = max_value /1000;

     %time % shit
    %value from trial 1
    tempStepVel = data_Noraxon(3).InsoleL.stepvelocity, ...
                  data_Noraxon(3).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(3).InsoleL.stepEndTime, ...
                  data_Noraxon(3).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(3).InsoleL.time(data_Noraxon(3).InsoleL.startIndx);
    stop = data_Noraxon(3).InsoleL.time(data_Noraxon(3).InsoleL.endIndx);
    length = stop - start;

    max_tPercent1 = ((max_t-start)/(length))*100;

    %value from trial 2
    tempStepVel = data_Noraxon(4).InsoleL.stepvelocity, ...
                  data_Noraxon(4).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(4).InsoleL.stepEndTime, ...
                  data_Noraxon(4).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(4).InsoleL.time(data_Noraxon(4).InsoleL.startIndx);
    stop = data_Noraxon(4).InsoleL.time(data_Noraxon(4).InsoleL.endIndx);
    length = stop - start;
    
    max_tPercent2 = ((max_t-start)/(length))*100;
            
    %median
    timeMaxVel = median([max_tPercent1,max_tPercent2]);

    nfs = struct('condition','nfs','median',median_value, 'iqr', iqr_value, 'max', max_value, 'timeatpeakvel',timeMaxVel);
%kss
    tempa = data_Noraxon(5).InsoleL.stepvelocity;
    tempb = data_Noraxon(5).InsoleR.stepvelocity;
    tempc = data_Noraxon(6).InsoleL.stepvelocity;
    tempd = data_Noraxon(6).InsoleR.stepvelocity;

    % Step 1: transpose
    tempa = tempa'; % Transpose to convert 1x16 to 16x1
    tempb = tempb'; % Transpose to convert 1x16 to 16x1
    tempc = tempc'; % Transpose to convert 1x16 to 16x1
    tempd = tempd'; % Transpose to convert 1x16 to 16x1
    
    % Concatenate vertically
    combined_array = [tempa; tempb; tempc; tempd]; 

    % Remove NaN values
    cleaned_array = rmmissing(combined_array);

    % Get unique values
    unique_array = unique(cleaned_array);

    combined_array = unique_array;

    % Step 4: Compute the calculations, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    step_velocities = [data_Noraxon(5).InsoleL.stepvelocity, ...
                  data_Noraxon(5).InsoleR.stepvelocity, ...
                  data_Noraxon(6).InsoleL.stepvelocity, ...
                  data_Noraxon(6).InsoleR.stepvelocity];

    step_endTime = [data_Noraxon(5).InsoleL.stepEndTime, ...
                  data_Noraxon(5).InsoleR.stepEndTime, ...
                  data_Noraxon(6).InsoleL.stepEndTime, ...
                  data_Noraxon(6).InsoleR.stepEndTime];

    %account for outliers
    q1 = median_value - iqr_value;
    q3 = median_value + iqr_value;

    iqr_ = q3 - q1;

    % Step 3: Define outliers as those below (Q1 - 1.5*IQR) or above (Q3 + 1.5*IQR)
    lower_bound = q1 - 1.5 * iqr_;
    upper_bound = q3 + 1.5 * iqr_;
    
    data  = step_velocities;
    % Step 4: Replace outliers with NaN
    data(data < lower_bound | data > upper_bound) = NaN;

    stepVelTime = [data; step_endTime];

    % Find the maximum value and its index
    [max_value, max_index] = max(stepVelTime(1,:));
    
    timeMaxVel = stepVelTime(2,max_index);

    median_value = median_value/1000; iqr_value = iqr_value/1000; max_value = max_value /1000;

     %time % shit
    %value from trial 1
    tempStepVel = data_Noraxon(5).InsoleL.stepvelocity, ...
                  data_Noraxon(5).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(5).InsoleL.stepEndTime, ...
                  data_Noraxon(5).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(5).InsoleL.time(data_Noraxon(5).InsoleL.startIndx);
    stop = data_Noraxon(5).InsoleL.time(data_Noraxon(5).InsoleL.endIndx);
    length = stop - start;

    max_tPercent1 = ((max_t-start)/(length))*100;

    %value from trial 2
    tempStepVel = data_Noraxon(6).InsoleL.stepvelocity, ...
                  data_Noraxon(6).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(6).InsoleL.stepEndTime, ...
                  data_Noraxon(6).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(6).InsoleL.time(data_Noraxon(6).InsoleL.startIndx);
    stop = data_Noraxon(6).InsoleL.time(data_Noraxon(6).InsoleL.endIndx);
    length = stop - start;
    
    max_tPercent2 = ((max_t-start)/(length))*100;
            
    %median
    timeMaxVel = median([max_tPercent1,max_tPercent2]);

    kss = struct('condition','kss','median',median_value, 'iqr', iqr_value, 'max', max_value, 'timeatpeakvel',timeMaxVel);
%kfs
    
    tempa = data_Noraxon(7).InsoleL.stepvelocity;
    tempb = data_Noraxon(7).InsoleR.stepvelocity;
    tempc = data_Noraxon(8).InsoleL.stepvelocity;
    tempd = data_Noraxon(8).InsoleR.stepvelocity;

    % Step 1: transpose
    tempa = tempa'; % Transpose to convert 1x16 to 16x1
    tempb = tempb'; % Transpose to convert 1x16 to 16x1
    tempc = tempc'; % Transpose to convert 1x16 to 16x1
    tempd = tempd'; % Transpose to convert 1x16 to 16x1
    
    % Concatenate vertically
    combined_array = [tempa; tempb; tempc; tempd]; 

    % Remove NaN values
    cleaned_array = rmmissing(combined_array);

    % Get unique values
    unique_array = unique(cleaned_array);

    combined_array = unique_array;

    % Step 4: Compute the calculations, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    step_velocities = [data_Noraxon(7).InsoleL.stepvelocity, ...
                  data_Noraxon(7).InsoleR.stepvelocity, ...
                  data_Noraxon(8).InsoleL.stepvelocity, ...
                  data_Noraxon(8).InsoleR.stepvelocity];

    step_endTime = [data_Noraxon(7).InsoleL.stepEndTime, ...
                  data_Noraxon(7).InsoleR.stepEndTime, ...
                  data_Noraxon(8).InsoleL.stepEndTime, ...
                  data_Noraxon(8).InsoleR.stepEndTime];

    %account for outliers
    q1 = median_value - iqr_value;
    q3 = median_value + iqr_value;

    iqr_ = q3 - q1;

    % Step 3: Define outliers as those below (Q1 - 1.5*IQR) or above (Q3 + 1.5*IQR)
    lower_bound = q1 - 1.5 * iqr_;
    upper_bound = q3 + 1.5 * iqr_;
    
    data  = step_velocities;
    % Step 4: Replace outliers with NaN
    data(data < lower_bound | data > upper_bound) = NaN;

    stepVelTime = [data; step_endTime];

    % Find the maximum value and its index
    [max_value, max_index] = max(stepVelTime(1,:));
    
    timeMaxVel = stepVelTime(2,max_index);

    median_value = median_value/1000; iqr_value = iqr_value/1000; max_value = max_value /1000;

     %time % shit
    %value from trial 1
    tempStepVel = data_Noraxon(7).InsoleL.stepvelocity, ...
                  data_Noraxon(7).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(7).InsoleL.stepEndTime, ...
                  data_Noraxon(7).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(7).InsoleL.time(data_Noraxon(7).InsoleL.startIndx);
    stop = data_Noraxon(7).InsoleL.time(data_Noraxon(7).InsoleL.endIndx);
    length = stop - start;

    max_tPercent1 = ((max_t-start)/(length))*100;

    %value from trial 2
    tempStepVel = data_Noraxon(8).InsoleL.stepvelocity, ...
                  data_Noraxon(8).InsoleR.stepvelocity;
    tempEndTime = data_Noraxon(8).InsoleL.stepEndTime, ...
                  data_Noraxon(8).InsoleR.stepEndTime;


    trial = [(tempStepVel);(tempEndTime)];
    
    [max_v, max_index] = max(trial(1,:));

    max_t = trial(2,max_index);

    start = data_Noraxon(8).InsoleL.time(data_Noraxon(8).InsoleL.startIndx);
    stop = data_Noraxon(8).InsoleL.time(data_Noraxon(8).InsoleL.endIndx);
    length = stop - start;
    
    max_tPercent2 = ((max_t-start)/(length))*100;
            
    %median
    timeMaxVel = median([max_tPercent1,max_tPercent2]);

    kfs = struct('condition','kfs','median',median_value, 'iqr', iqr_value, 'max', max_value, 'timeatpeakvel',timeMaxVel);

%
data_velocity = [nss, nfs, kss, kfs];
end