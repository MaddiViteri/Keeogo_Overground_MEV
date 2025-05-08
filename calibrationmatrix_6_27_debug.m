function [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha] = calibrationmatrix_6_27_debug(static_markerStruct,markerStruct)
    % Save the input structures for inspection
    save('debug_inputs.mat', 'static_markerStruct', 'markerStruct');

    % List of required markers
    requiredMarkers = {'LLatMal', 'RLatMal', 'LMedMal', 'RMedMal'};  % Add others as needed

    % Check and allow manual input for any missing marker
    for i = 1:length(requiredMarkers)
        fieldName = requiredMarkers{i};
        if ~isfield(static_markerStruct, fieldName)
            fprintf('\nMarker "%s" is missing.\n', fieldName);
            userInput = input(sprintf('Enter field name in static_markerStruct to use for "%s":\n', fieldName), 's');
                try
                    static_markerStruct.(fieldName) = static_markerStruct.(userInput);
                catch ME
                    error('Failed to assign marker "%s": %s', fieldName, ME.message);
                end
        end
    end

    % Call original function with fixed inputs
    [right_globe2ana,left_globe2ana,right_seg_omega,left_seg_omega,right_seg_alpha,left_seg_alpha,right_ang,right_omega,right_alpha,left_alpha] = calibrationmatrix_6_27(static_markerStruct,markerStruct);
end