%% MemToolbox fitting: participant x temporal location
% Model: Orientation(WithBias(StandardMixtureModel()), [1 3])
% Transformation: processedError = asind(rawError / 80)
% Estimates: maximum a posteriori estimates (fitResult.maxPosterior)
% Output: SubjectID, td, mu, sd, g
% Note: mu and sd are in degrees on the transformed scale
% Exclude participants 12 and 13 due to non-convergent fits

clear; clc;
rng(20260314, 'twister');   % Fixed seed for reproducibility

%% Paths
% Put this script next to a folder named "signederror_eachsub"
scriptPath = mfilename('fullpath');
scriptDir  = fileparts(scriptPath);

% When the script is run section-by-section, mfilename may be empty.
if isempty(scriptDir)
    scriptDir = pwd;
end

dataPath = fullfile(scriptDir, 'signederror_eachsub');
savePath = fullfile(scriptDir, 'memresults.xlsx');

% Add MemToolbox to the MATLAB path first if necessary:
% addpath(genpath('E:\path\to\MemToolbox'));

assert(exist('MemFit', 'file') == 2, ...
    'MemToolbox is not on the MATLAB path.');
assert(isfolder(dataPath), ...
    'Data folder does not exist: %s', dataPath);

%% Analysis settings
tdLevels    = [0.2, 0.4, 0.6, 0.8];
tdTolerance = 1e-8;

model      = Orientation(WithBias(StandardMixtureModel()), [1 3]);
paramNames = model.paramNames;

muIdx = find(strcmpi(paramNames, 'mu'), 1);
sdIdx = find(strcmpi(paramNames, 'sd'), 1);
gIdx  = find(strcmpi(paramNames, 'g'),  1);

assert(~isempty(muIdx) && ~isempty(sdIdx) && ~isempty(gIdx), ...
    'Could not identify mu, sd, and g in model.paramNames.');

%% Input files
fileList = dir(fullfile(dataPath, '*.xlsx'));

if isempty(fileList)
    error('No .xlsx files found in: %s', dataPath);
end

% Fixed file order helps reproducibility.
[~, order] = sort(lower({fileList.name}));
fileList = fileList(order);

%% Results
results = {};

%% Fit each participant and td level
for f = 1:numel(fileList)

    fileName = fileList(f).name;
    filePath = fullfile(dataPath, fileName);
    [~, subID, ~] = fileparts(fileName);

    % Excluded participants
    if startsWith(fileName, '12_') || startsWith(fileName, '13_')
        fprintf('\nSkipped excluded participant: %s\n', fileName);
        continue;
    end

    fprintf('\n========================================\n');
    fprintf('Participant: %s\n', subID);
    fprintf('========================================\n');

    %% Read data
    try
        data = readmatrix(filePath);
    catch ME
        warning('Could not read %s: %s', fileName, ME.message);
        continue;
    end

    if size(data, 2) < 8
        warning('%s contains fewer than 8 columns and was skipped.', fileName);
        continue;
    end

    tdCol    = data(:, 7);
    errorCol = data(:, 8);

    %% Fit each td level
    for t = 1:numel(tdLevels)

        td = tdLevels(t);

        % Avoid direct equality comparisons for floating-point values.
        conditionMask = isfinite(tdCol) & abs(tdCol - td) < tdTolerance;

        if ~any(conditionMask)
            warning('Participant %s has no data at td=%.1f.', subID, td);
            continue;
        end

        rawError = errorCol(conditionMask);

        % Remove missing or non-finite values.
        rawError = rawError(isfinite(rawError));

        % Values outside this range make asind(error/80) invalid.
        if any(abs(rawError) > 80 + 1e-10)
            warning(['Participant %s, td=%.1f contains error values ', ...
                     'outside [-80, 80] and was skipped.'], subID, td);
            continue;
        end

        processedError = asind(rawError / 80);
        processedError = processedError(isfinite(processedError));
        
        memData = struct();
        memData.errors = processedError(:)';

        %% MemToolbox fit
        try
            fitResult = MemFit(memData, model, 'Verbosity', 0);
            mapParams = fitResult.maxPosterior;

            if any(~isfinite(mapParams))
                warning('Participant %s, td=%.1f returned NaN or Inf.', ...
                    subID, td);
                continue;
            end

            muVal = mapParams(muIdx);
            sdVal = mapParams(sdIdx);
            gVal  = mapParams(gIdx);

            fprintf('  td=%.1f | mu=%.6f | sd=%.6f | g=%.6f\n', ...
                td, muVal, sdVal, gVal);

            results(end+1, :) = {subID, td, muVal, sdVal, gVal}; %#ok<SAGROW>

        catch ME
            warning('MemFit failed for participant %s, td=%.1f: %s', ...
                subID, td, ME.message);
        end
    end
end

%% Save results
if isempty(results)
    warning('No fits were completed successfully.');
else
    resultTable = cell2table(results, ...
        'VariableNames', {'SubjectID', 'td', 'mu', 'sd', 'g'});

    resultTable = sortrows(resultTable, {'SubjectID', 'td'});
    writetable(resultTable, savePath);

    fprintf('\nSuccessful fits: %d\n', height(resultTable));
    disp(resultTable);
    fprintf('Results saved to:\n%s\n', savePath);
end
