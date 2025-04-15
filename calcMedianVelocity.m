function median_velocity = calcMedianVelocity(timeArrays, velocityArrays)
   % Determine a common time vector based on the full range of time data
    minTime = min(cellfun(@min, timeArrays)); % Earliest start time
    maxTime = max(cellfun(@max, timeArrays)); % Latest end time
    commonTime = linspace(minTime, maxTime, 100); % 100 evenly spaced points
    
    % Interpolate velocities onto the common time grid
    numTrials = length(timeArrays);
    interpolatedVelocities = zeros(numTrials, length(commonTime));
    
    for i = 1:numTrials
        interpolatedVelocities(i, :) = interp1(timeArrays{i}, velocityArrays{i}, commonTime, 'linear', 'extrap');
    end
    
    % Compute the overall median velocity across all trials and time points
    overall_median = median(interpolatedVelocities(:),"omitmissing");
    
    median_velocity = overall_median;
end