function [data_Noraxon] = calcVel (data_Noraxon)

    for m = 1:8
        % Left Insole
        startIndxL = data_Noraxon(m).InsoleL.startIndx;
        endIndxL = data_Noraxon(m).InsoleL.endIndx;
        
        % Extract COP_X and COP_Y data for the relevant time range (Left Insole)
        copY_L = data_Noraxon(m).InsoleL.COP_Y(startIndxL:endIndxL);
        copX_L = data_Noraxon(m).InsoleL.COP_X(startIndxL:endIndxL);

        %convert from mm to m
        copY_L = copY_L./1000;
        copX_L = copX_L./1000;
        
        % Calculate velocity components (using diff to approximate)
        velX_L = diff(copX_L);  % Velocity in X direction (using difference)
        velY_L = diff(copY_L);  % Velocity in Y direction (using difference)
        
        % Combine velocity components to calculate the total velocity
        vel_L = sqrt(velX_L.^2 + velY_L.^2);  % Total velocity (magnitude of the velocity vector)
        data_Noraxon(m).InsoleL.vel = vel_L;
        
        % Omit NaNs in the velocity data
        valid_vel_L = vel_L(~isnan(vel_L));  % Remove NaNs from velocity data
        
        % Calculate statistics for Left Insole
        max_vel_L = max(valid_vel_L);
        median_vel_L = median(valid_vel_L);
        iqr_vel_L = iqr(valid_vel_L);
        
        % Time at maximum velocity
        [~, max_vel_idx_L] = max(valid_vel_L);
        time_at_max_vel_L = data_Noraxon(m).InsoleL.time(startIndxL + max_vel_idx_L);  % Corresponding time index for max velocity
        
        % Store the results in the data structure
        data_Noraxon(m).InsoleL.max_velocity = max_vel_L;
        data_Noraxon(m).InsoleL.median_velocity = median_vel_L;
        data_Noraxon(m).InsoleL.iqr_velocity = iqr_vel_L;
        data_Noraxon(m).InsoleL.time_at_max_velocity = time_at_max_vel_L;
        
        
        % Right Insole
        startIndxR = data_Noraxon(m).InsoleR.startIndx;
        endIndxR = data_Noraxon(m).InsoleR.endIndx;
        
        % Extract COP_X and COP_Y data for the relevant time range (Right Insole)
        copY_R = data_Noraxon(m).InsoleR.COP_Y(startIndxR:endIndxR);
        copX_R = data_Noraxon(m).InsoleR.COP_X(startIndxR:endIndxR);

        %convert from mm to m
        copY_R = copY_R./1000;
        copX_R = copX_R./1000;
        
        % Calculate velocity components (using diff to approximate)
        velX_R = diff(copX_R);  % Velocity in X direction (using difference)
        velY_R = diff(copY_R);  % Velocity in Y direction (using difference)
        
        % Combine velocity components to calculate the total velocity
        vel_R = sqrt(velX_R.^2 + velY_R.^2);  % Total velocity (magnitude of the velocity vector)
        data_Noraxon(m).InsoleR.vel = vel_R;
        
        % Omit NaNs in the velocity data
        valid_vel_R = vel_R(~isnan(vel_R));  % Remove NaNs from velocity data
        
        % Calculate statistics for Right Insole
        max_vel_R = max(valid_vel_R);
        median_vel_R = median(valid_vel_R);
        iqr_vel_R = iqr(valid_vel_R);
        
        % Time at maximum velocity
        [~, max_vel_idx_R] = max(valid_vel_R);
        time_at_max_vel_R = data_Noraxon(m).InsoleR.time(startIndxR + max_vel_idx_R);  % Corresponding time index for max velocity
        
        % Store the results in the data structure
        data_Noraxon(m).InsoleR.max_velocity = max_vel_R;
        data_Noraxon(m).InsoleR.median_velocity = median_vel_R;
        data_Noraxon(m).InsoleR.iqr_velocity = iqr_vel_R;
        data_Noraxon(m).InsoleR.time_at_max_velocity = time_at_max_vel_R;
        
    end % End of for loop

end % End of function
