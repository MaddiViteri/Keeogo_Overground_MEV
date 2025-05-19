% compare_txt_duplicates_mac_struct_check.m
% Reads .txt files into a struct, detects duplicates, renames them,
% and verifies PlatformRotation against filename label.

% Select folder
folderPath = uigetdir(pwd, 'Select folder with .txt files');
if folderPath == 0
    error('No folder selected.');
end

% Define expected PlatformRotation values
expectedRotation.downhill = -8.6;
expectedRotation.neutral  =  0.0;
expectedRotation.uphill   =  8.6;

% Get list of .txt files
txtFiles = dir(fullfile(folderPath, '*.txt'));
nFiles = length(txtFiles);

% Containers
fileHashes = containers.Map();     % MD5 hash -> filename
duplicates = {};                   % Duplicate file log
fileStruct = struct();             % Structure of tables
rotationMismatches = {};           % Store mismatches

% Loop through files
for i = 1:nFiles
    filename = txtFiles(i).name;
    filePath = fullfile(folderPath, filename);
    
    % Create a valid struct field name
    [~, nameOnly, ~] = fileparts(filename);
    structField = matlab.lang.makeValidName(nameOnly);

    try
        % Read file as table
        tbl = readtable(filePath, 'PreserveVariableNames', true);
        fileStruct.(structField) = tbl;
    catch
        warning('⚠️ Could not read %s as table. Skipping.', filename);
        continue;
    end

    % Check PlatformRotation
    platformValue = NaN;
    if any(strcmp('PlatformRotation', tbl.Properties.VariableNames))
        platformValue = unique(tbl.PlatformRotation);
        platformValue = platformValue(1); % Just take first if multiple
    else
        warning('⚠️ File "%s" missing PlatformRotation column.', filename);
        continue;
    end

    % Determine expected rotation based on filename
    matched = false;
    for condition = fieldnames(expectedRotation)'
        condStr = condition{1};
        if contains(lower(filename), condStr)
            matched = true;
            expectedVal = expectedRotation.(condStr);
            if abs(platformValue - expectedVal) > 0.01
                fprintf('❌ Mismatch: %s has PlatformRotation = %.1f (expected %.1f)\n', ...
                        filename, platformValue, expectedVal);
                rotationMismatches{end+1, 1} = filename; %#ok<AGROW>
                rotationMismatches{end, 2} = platformValue;
                rotationMismatches{end, 3} = expectedVal;
            end
            break;
        end
    end

    if ~matched
        warning('⚠️ Could not determine condition from filename: %s', filename);
    end

    % Check for duplicates
    tableString = evalc('disp(tbl)');
    hash = getMD5(tableString);
    if isKey(fileHashes, hash)
        original = fileHashes(hash);
        fprintf('🔁 Duplicate found:\n  %s\n  --> is a copy of: %s\n\n', ...
                filename, original);
        duplicates{end+1, 1} = filename; %#ok<AGROW>
        duplicates{end, 2} = original;

        % Rename file
        newName = ['DUPLICATE_' filename];
        newPath = fullfile(folderPath, newName);
        counter = 1;
        while exist(newPath, 'file')
            newName = ['DUPLICATE_' num2str(counter) '_' filename];
            newPath = fullfile(folderPath, newName);
            counter = counter + 1;
        end
        movefile(filePath, newPath);
    else
        fileHashes(hash) = filename;
    end
end

% Final reports
if isempty(duplicates)
    disp('✅ No duplicate files found.');
else
    disp('📁 Duplicate files renamed:');
    disp(duplicates);
end

if isempty(rotationMismatches)
    disp('✅ All PlatformRotation values match their filenames.');
else
    fprintf('\n🛑 PlatformRotation mismatches found:\n');
    for i = 1:size(rotationMismatches, 1)
        fprintf('  %s → %.1f (expected %.1f)\n', ...
                rotationMismatches{i,1}, rotationMismatches{i,2}, rotationMismatches{i,3});
    end
end

disp('🗂️ All tables stored in: fileStruct');
