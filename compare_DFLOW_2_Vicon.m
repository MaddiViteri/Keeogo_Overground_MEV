% File: check_condition_files.m

% Choose folder
folderPath = uigetdir(pwd, 'Select folder with .txt files');
if folderPath == 0
    error('No folder selected.');
end

% Define conditions and expected base names
conditions = {'downhill', 'neutral', 'uphill'};
expectedSuffixes = {'static', 'control', 'control01', 'keeogo', 'keeogo01'};

% Get all .txt files in the folder
fileList = dir(fullfile(folderPath, '*.txt'));
fileNames = lower({fileList.name});  % Lowercase for case-insensitive matching
fileNames = string(fileNames);

% Build expected full filenames
expectedFiles = strings(0);
for i = 1:length(conditions)
    for j = 1:length(expectedSuffixes)
        expectedFiles(end+1) = lower(sprintf('%s_%s.txt', conditions{i}, expectedSuffixes{j}));
    end
end

% Check which expected files are present
fprintf('\n📂 Checking files in: %s\n\n', folderPath);
foundFiles = [];
missingFiles = [];

for i = 1:length(expectedFiles)
    fname = expectedFiles(i);
    if any(fileNames == fname)
        fprintf('  ✅ Found: %s\n', fname);
        foundFiles(end+1) = fname; %#ok<AGROW>
    else
        fprintf('  ❌ Missing: %s\n', fname);
        missingFiles(end+1) = fname; %#ok<AGROW>
    end
end

% Check for unexpected files
unexpectedFiles = setdiff(fileNames, expectedFiles);

if ~isempty(unexpectedFiles)
    fprintf('\n⚠️ Unexpected files found (not part of expected naming scheme):\n');
    disp(unexpectedFiles');
else
    fprintf('\n✅ No unexpected files found.\n');
end

% Summary
fprintf('\n📊 Summary:\n');
fprintf('  Expected files: %d\n', length(expectedFiles));
fprintf('  Found: %d\n', length(foundFiles));
fprintf('  Missing: %d\n', length(missingFiles));
fprintf('  Unexpected: %d\n', length(unexpectedFiles));
