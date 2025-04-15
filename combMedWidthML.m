function [nss_ML, nfs_ML, kss_ML, kfs_ML] = combMedWidthML (data_Noraxon)
    tempa = data_Noraxon(1).InsoleL.width_ML;
    tempb = data_Noraxon(1).InsoleR.width_ML;
    tempc = data_Noraxon(2).InsoleL.width_ML;
    tempd = data_Noraxon(2).InsoleR.width_ML;
    
    % Step 1: Find the maximum length
    maxLength = max([length(tempa), length(tempb), length(tempc), length(tempd)]);
    
    % Step 2: Pad the arrays with NaNs (or zeros) to match the maximum length
    tempa_padded = [tempa, nan(1, maxLength - length(tempa))];
    tempb_padded = [tempb, nan(1, maxLength - length(tempb))];
    tempc_padded = [tempc, nan(1, maxLength - length(tempc))];
    tempd_padded = [tempd, nan(1, maxLength - length(tempd))];
    
    % Step 3: Concatenate the padded arrays into one combined array
    combined_array = [tempa_padded, tempb_padded, tempc_padded, tempd_padded];
    combined_array = rmmissing(combined_array);

    % Step 4: Compute the median, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    nss_ML = struct('median',median_value,'IQR',iqr_value);
    
    tempa = data_Noraxon(3).InsoleL.width_ML;
    tempb = data_Noraxon(3).InsoleR.width_ML;
    tempc = data_Noraxon(4).InsoleL.width_ML;
    tempd = data_Noraxon(4).InsoleR.width_ML;
    
    % Step 1: Find the maximum length
    maxLength = max([length(tempa), length(tempb), length(tempc), length(tempd)]);
    
    % Step 2: Pad the arrays with NaNs (or zeros) to match the maximum length
    tempa_padded = [tempa, nan(1, maxLength - length(tempa))];
    tempb_padded = [tempb, nan(1, maxLength - length(tempb))];
    tempc_padded = [tempc, nan(1, maxLength - length(tempc))];
    tempd_padded = [tempd, nan(1, maxLength - length(tempd))];
    
    % Step 3: Concatenate the padded arrays into one combined array
    combined_array = [tempa_padded, tempb_padded, tempc_padded, tempd_padded];
    combined_array = rmmissing(combined_array);

    % Step 4: Compute the median, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    nfs_ML = struct('median',median_value,'IQR',iqr_value);
    
    tempa = data_Noraxon(5).InsoleL.width_ML;
    tempb = data_Noraxon(5).InsoleR.width_ML;
    tempc = data_Noraxon(6).InsoleL.width_ML;
    tempd = data_Noraxon(6).InsoleR.width_ML;
    
    % Step 1: Find the maximum length
    maxLength = max([length(tempa), length(tempb), length(tempc), length(tempd)]);
    
    % Step 2: Pad the arrays with NaNs (or zeros) to match the maximum length
    tempa_padded = [tempa, nan(1, maxLength - length(tempa))];
    tempb_padded = [tempb, nan(1, maxLength - length(tempb))];
    tempc_padded = [tempc, nan(1, maxLength - length(tempc))];
    tempd_padded = [tempd, nan(1, maxLength - length(tempd))];
    
    % Step 3: Concatenate the padded arrays into one combined array
    combined_array = [tempa_padded, tempb_padded, tempc_padded, tempd_padded];
    combined_array = rmmissing(combined_array);

    % Step 4: Compute the median, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    kss_ML = struct('median',median_value,'IQR',iqr_value);
    
    tempa = data_Noraxon(7).InsoleL.width_ML;
    tempb = data_Noraxon(7).InsoleR.width_ML;
    tempc = data_Noraxon(8).InsoleL.width_ML;
    tempd = data_Noraxon(8).InsoleR.width_ML;
    
    % Step 1: Find the maximum length
    maxLength = max([length(tempa), length(tempb), length(tempc), length(tempd)]);
    
    % Step 2: Pad the arrays with NaNs (or zeros) to match the maximum length
    tempa_padded = [tempa, nan(1, maxLength - length(tempa))];
    tempb_padded = [tempb, nan(1, maxLength - length(tempb))];
    tempc_padded = [tempc, nan(1, maxLength - length(tempc))];
    tempd_padded = [tempd, nan(1, maxLength - length(tempd))];
    
    % Step 3: Concatenate the padded arrays into one combined array
    combined_array = [tempa_padded, tempb_padded, tempc_padded, tempd_padded];
    combined_array = rmmissing(combined_array);

    % Step 4: Compute the median, omitting NaNs
    median_value = median(combined_array, 'omitnan');
    iqr_value = iqr(combined_array);
    
    kfs_ML = struct('median',median_value,'IQR',iqr_value);
end