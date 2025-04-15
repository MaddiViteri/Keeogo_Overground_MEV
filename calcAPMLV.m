function [data_Noraxon] = calcAPMLV(data_Noraxon)

    for m = 1:8
        
        % Left Insole
        startIndxL = data_Noraxon(m).InsoleL.startIndx;
        endIndxL = data_Noraxon(m).InsoleL.endIndx;
        
        % Extract COP and time data for Left Insole
        copY_L = data_Noraxon(m).InsoleL.COP_Y(startIndxL:endIndxL);
        copX_L = data_Noraxon(m).InsoleL.COP_X(startIndxL:endIndxL);
        time_L = data_Noraxon(m).InsoleL.time(startIndxL:endIndxL); % Assuming time is available
        
        % Initialize arrays
        pairedLengths_L = [];
        pairedWidths_L = [];
        stepVelocities_L = [];
        stependtime_L = [];
        
        % Find valid step segments using non-NaN values
        nan_indicesL = isnan(copY_L);
        valid_segmentsL = find(~nan_indicesL);
        nan_indicesW = isnan(copX_L);
        valid_segmentsW = find(~nan_indicesW);
        
        % Extract contiguous step segments
        segmentsL = {};
        segmentsW = {};
        current_segmentL = [];
        current_segmentW = [];
        
        for i = 1:length(valid_segmentsL)
            if isempty(current_segmentL) || valid_segmentsL(i) == valid_segmentsL(i-1) + 1
                current_segmentL(end+1) = valid_segmentsL(i);
            else
                segmentsL{end+1} = current_segmentL;
                current_segmentL = valid_segmentsL(i);
            end
        end
        if ~isempty(current_segmentL)
            segmentsL{end+1} = current_segmentL;
        end
        
        for i = 1:length(valid_segmentsW)
            if isempty(current_segmentW) || valid_segmentsW(i) == valid_segmentsW(i-1) + 1
                current_segmentW(end+1) = valid_segmentsW(i);
            else
                segmentsW{end+1} = current_segmentW;
                current_segmentW = valid_segmentsW(i);
            end
        end
        if ~isempty(current_segmentW)
            segmentsW{end+1} = current_segmentW;
        end
        
        % Process each valid step segment
        for i = 1:length(segmentsL)
            valid_data_Y = copY_L(segmentsL{i});
            valid_data_X = copX_L(segmentsW{i});
            valid_time = time_L(segmentsL{i});
            
            if ~isempty(valid_data_Y) && ~isempty(valid_data_X) && numel(valid_time) > 1
                % Compute AP & ML lengths
                length_val = max(valid_data_Y) - min(valid_data_Y);
                width_val = max(valid_data_X) - min(valid_data_X);
                
                if length_val > 0, pairedLengths_L(end+1) = length_val; end
                if width_val > 0, pairedWidths_L(end+1) = width_val; end
                
                % Compute step displacement (Euclidean distance)
                dx = valid_data_X(end) - valid_data_X(1);
                dy = valid_data_Y(end) - valid_data_Y(1);
                displacement = sqrt(dx^2 + dy^2);
                
                % Compute step time
                step_time = valid_time(end) - valid_time(1);
                stependtime_L (end+1) = valid_time(end);
                
                % Compute velocity (Displacement / Time)
                if step_time > 0
                    stepVelocities_L(end+1) = displacement / step_time;
                end
            end
        end
        
        % Store results in dataset
        data_Noraxon(m).InsoleL.length_AP = pairedLengths_L;
        data_Noraxon(m).InsoleL.width_ML = pairedWidths_L;
        data_Noraxon(m).InsoleL.stepvelocity = stepVelocities_L;
        data_Noraxon(m).InsoleL.stepEndTime = stependtime_L;

        % Compute median values
        data_Noraxon(m).InsoleL.mAP = median(pairedLengths_L);
        data_Noraxon(m).InsoleL.mML = median(pairedWidths_L);
        data_Noraxon(m).InsoleL.mVel = median(stepVelocities_L);
        
        % Repeat for Right Insole
        startIndxR = data_Noraxon(m).InsoleR.startIndx;
        endIndxR = data_Noraxon(m).InsoleR.endIndx;
        
        copY_R = data_Noraxon(m).InsoleR.COP_Y(startIndxR:endIndxR);
        copX_R = data_Noraxon(m).InsoleR.COP_X(startIndxR:endIndxR);
        time_R = data_Noraxon(m).InsoleR.time(startIndxR:endIndxR);
        
        pairedLengths_R = [];
        pairedWidths_R = [];
        stepVelocities_R = [];
        stependtime_R =[];
        
        nan_indicesL_R = isnan(copY_R);
        valid_segmentsL_R = find(~nan_indicesL_R);
        nan_indicesW_R = isnan(copX_R);
        valid_segmentsW_R = find(~nan_indicesW_R);
        
        segmentsL_R = {};
        segmentsW_R = {};
        current_segmentL_R = [];
        current_segmentW_R = [];
        
        for i = 1:length(valid_segmentsL_R)
            if isempty(current_segmentL_R) || valid_segmentsL_R(i) == valid_segmentsL_R(i-1) + 1
                current_segmentL_R(end+1) = valid_segmentsL_R(i);
            else
                segmentsL_R{end+1} = current_segmentL_R;
                current_segmentL_R = valid_segmentsL_R(i);
            end
        end
        if ~isempty(current_segmentL_R)
            segmentsL_R{end+1} = current_segmentL_R;
        end
        
        for i = 1:length(valid_segmentsW_R)
            if isempty(current_segmentW_R) || valid_segmentsW_R(i) == valid_segmentsW_R(i-1) + 1
                current_segmentW_R(end+1) = valid_segmentsW_R(i);
            else
                segmentsW_R{end+1} = current_segmentW_R;
                current_segmentW_R = valid_segmentsW_R(i);
            end
        end
        if ~isempty(current_segmentW_R)
            segmentsW_R{end+1} = current_segmentW_R;
        end
        
        for i = 1:length(segmentsL_R)
            valid_data_Y = copY_R(segmentsL_R{i});
            valid_data_X = copX_R(segmentsW_R{i});
            valid_time = time_R(segmentsL_R{i});
            
            if ~isempty(valid_data_Y) && ~isempty(valid_data_X) && numel(valid_time) > 1
                length_val = max(valid_data_Y) - min(valid_data_Y);
                width_val = max(valid_data_X) - min(valid_data_X);
                
                if length_val > 0, pairedLengths_R(end+1) = length_val; end
                if width_val > 0, pairedWidths_R(end+1) = width_val; end
                
                dx = valid_data_X(end) - valid_data_X(1);
                dy = valid_data_Y(end) - valid_data_Y(1);
                displacement = sqrt(dx^2 + dy^2);
                
                step_time = valid_time(end) - valid_time(1);
                stependtime_R (end+1) = valid_time(end);
                
                if step_time > 0
                    stepVelocities_R(end+1) = displacement / step_time;
                end
            end
        end
        
        data_Noraxon(m).InsoleR.length_AP = pairedLengths_R;
        data_Noraxon(m).InsoleR.width_ML = pairedWidths_R;
        data_Noraxon(m).InsoleR.stepvelocity = stepVelocities_R;
        data_Noraxon(m).InsoleR.stepEndTime = stependtime_R;


        data_Noraxon(m).InsoleR.mAP = median(pairedLengths_R);
        data_Noraxon(m).InsoleR.mML = median(pairedWidths_R);
        data_Noraxon(m).InsoleR.mVel = median(stepVelocities_R);
    end

end
