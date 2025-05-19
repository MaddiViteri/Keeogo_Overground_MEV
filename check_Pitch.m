% compare_txt_duplicates_mac_struct_check_rename.m
% Reads .txt files into a struct, detects duplicates, renames them,
% verifies PlatformRotation vs filename label,
% renames mismatched or invalid files with appropriate prefixes.

% Select folder
folderPath = uigetdir(pwd, 'Select folder with .txt files');
if folderPath == 0
    error('No folder selected.');
end

disp(folderPath);

% Expected PlatformRotation values
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
    
    % Make valid struct field name
    [~, nameOnly, ext] = fileparts(filename);
    structField = matlab.lang.makeValidName(nameOnly);

    try
        % Read as table
        tbl = readtable(filePath, 'PreserveVariableNames', true);
        fileStruct.(structField) = tbl;
    catch
        warning('⚠️ Could not read %s as table. Skipping.', filename);
        continue;
    end

    % Check PlatformRotation column presence
    if ~any(strcmp('PlatformRotation', tbl.Properties.VariableNames))
        warning('⚠️ File "%s" missing PlatformRotation column.', filename);

        % Rename file with INVALID_ prefix
        newName = ['INVALID_' filename];
        newPath = fullfile(folderPath, newName);
        counter = 1;
        while exist(newPath, 'file')
            newName = sprintf('INVALID_%d_%s', counter, filename);
            newPath = fullfile(folderPath, newName);
            counter = counter + 1;
        end
        movefile(filePath, newPath);
        continue; % Skip further processing for this file
    end

    % Get PlatformRotation value (first unique value)
    platformValue = unique(tbl.PlatformRotation);
    platformValue = platformValue(1);

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

                % Rename mismatched file with prefix ROT_<value>_
                prefix = sprintf('ROT_%.1f_', platformValue);
                newName = [prefix filename];
                newPath = fullfile(folderPath, newName);
                counter = 1;
                while exist(newPath, 'file')
                    newName = sprintf('%s%d_%s', prefix, counter, filename);
                    newPath = fullfile(folderPath, newName);
                    counter = counter + 1;
                end
                movefile(filePath, newPath);
            end
            break;
        end
    end

    if ~matched
        warning('⚠️ Could not determine condition from filename: %s', filename);
    end

    % Compute hash on current file content (after possible renaming)
    % Note: If renamed, update filename & filePath accordingly
    if exist(newPath, 'file')
        % File was renamed
        [~, filename] = fileparts(newName);
        filename = [filename ext];
        filePath = newPath;
    end

    % Reload table after rename (optional, to keep consistency)
    try
        tbl = readtable(filePath, 'PreserveVariableNames', true);
        fileStruct.(matlab.lang.makeValidName(nameOnly)) = tbl;
    catch
        warning('⚠️ Could not reload renamed file %s.', filename);
        % We'll proceed anyway
    end

    % Compute MD5 hash string of the table string
    tableString = evalc('disp(tbl)');
    hash = getMD5(tableString);

    if isKey(fileHashes, hash)
        original = fileHashes(hash);
        fprintf('🔁 Duplicate found:\n  %s\n  --> is a copy of: %s\n\n', ...
                filename, original);
        duplicates{end+1, 1} = filename; %#ok<AGROW>
        duplicates{end, 2} = original;

        % Rename duplicate file with DUPLICATE_ prefix
        newNameDup = ['DUPLICATE_' filename];
        newPathDup = fullfile(folderPath, newNameDup);
        counter = 1;
        while exist(newPathDup, 'file')
            newNameDup = sprintf('DUPLICATE_%d_%s', counter, filename);
            newPathDup = fullfile(folderPath, newNameDup);
            counter = counter + 1;
        end
        movefile(filePath, newPathDup);
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

% Create summary table of filenames and PlatformRotation values
summaryData = {};
summaryHeader = {'Filename', 'PlatformRotation'};

allFields = fieldnames(fileStruct);
for i = 1:numel(allFields)
    field = allFields{i};
    tbl = fileStruct.(field);
    if any(strcmp('PlatformRotation', tbl.Properties.VariableNames))
        platformVal = unique(tbl.PlatformRotation);
        summaryData{end+1, 1} = field; %#ok<AGROW>
        summaryData{end, 2} = platformVal(1);
    else
        summaryData{end+1, 1} = field; %#ok<AGROW>
        summaryData{end, 2} = NaN; % Indicate missing column
    end
end

summaryTable = cell2table(summaryData, 'VariableNames', summaryHeader);
disp('📊 Summary of PlatformRotation values:');
disp(summaryTable);

disp('🗂️ All tables stored in: fileStruct');
