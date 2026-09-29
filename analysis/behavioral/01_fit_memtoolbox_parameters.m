%% Fit MemToolbox parameters for the three behavioral groups
%
% Required software:
%   - MATLAB
%   - MemToolbox (Suchow et al., 2013) available on the MATLAB path
%
% Model:
%   Orientation(WithBias(StandardMixtureModel()), [1 3])
%
% The model is fitted separately for every participant and temporal
% location (20%, 40%, 60%, and 80%). Signed errors in percentage points
% are transformed with asind(error / 80) before fitting.
%
% This script intentionally does not create convergence-status variables,
% boundary flags, or convergence log files. It writes fitted parameters to
% results/behavioral/memtoolbox_parameters_reproduced.csv and never
% overwrites the cleaned source data.

clear;
clc;

%% Reproducibility and paths
rng(20260917, 'twister');

script_path = mfilename('fullpath');
if isempty(script_path)
    error(['Unable to determine the script location. Run this file as a ', ...
        'saved MATLAB script.']);
end

script_dir = fileparts(script_path);
repo_root = fileparts(fileparts(script_dir));
output_dir = fullfile(repo_root, 'results', 'behavioral');

if ~isfolder(output_dir)
    mkdir(output_dir);
end

if exist('MemFit', 'file') ~= 2 || ...
        exist('Orientation', 'file') ~= 2 || ...
        exist('WithBias', 'file') ~= 2 || ...
        exist('StandardMixtureModel', 'file') ~= 2
    error(['MemToolbox was not found on the MATLAB path. Download and ', ...
        'add MemToolbox to the path before running this script.']);
end

%% Input files and manuscript-defined parameter-analysis sample
group_labels = {'fMRI', 'replication', 'control'};

input_files = {
    fullfile(repo_root, 'data', 'behavioral', 'fMRI_group', ...
        'fmri_group_trials.csv')
    fullfile(repo_root, 'data', 'behavioral', 'replication_group', ...
        'replication_group_trials.csv')
    fullfile(repo_root, 'data', 'behavioral', 'control_group', ...
        'control_group_trials.csv')
};

location_columns = {
    'temporal_location(%)'
    'target_location'
    'temporal_location(%)'
};

% The manuscript reports MemToolbox parameter analyses for 16 of the 18
% fMRI participants. These two documented exclusions are applied directly
% so that the analysis sample is explicit and reproducible.
excluded_subjects = {
    {'sub-12', 'sub-13'}
    {}
    {}
};

temporal_locations = [20, 40, 60, 80];
minimum_trials_per_fit = 5;

%% Define the MemToolbox model and parameter positions
model = Orientation(WithBias(StandardMixtureModel()), [1 3]);
parameter_names = lower(string(model.paramNames));

mu_index = find(parameter_names == "mu", 1);
sigma_index = find(parameter_names == "sd" | ...
    parameter_names == "sigma", 1);
g_index = find(parameter_names == "g", 1);

if isempty(mu_index) || isempty(sigma_index) || isempty(g_index)
    error('Could not identify mu, sigma/sd, and g in model.paramNames.');
end

%% Fit each participant and temporal location
fit_rows = cell(0, 6);

for group_index = 1:numel(group_labels)
    group_label = group_labels{group_index};
    input_file = input_files{group_index};
    location_column = location_columns{group_index};

    if ~isfile(input_file)
        error('Required input file was not found: %s', input_file);
    end

    import_options = detectImportOptions(input_file, ...
        'VariableNamingRule', 'preserve');
    trial_data = readtable(input_file, import_options);

    required_columns = {'subject', location_column, 'signed_error'};
    missing_columns = setdiff(required_columns, ...
        trial_data.Properties.VariableNames);
    if ~isempty(missing_columns)
        error('Input file %s is missing required columns: %s', ...
            input_file, strjoin(missing_columns, ', '));
    end

    subject_values = string(trial_data.('subject'));
    location_values = double(trial_data.(location_column));
    error_values = double(trial_data.('signed_error'));

    available_locations = sort(unique(location_values(isfinite(location_values))))';
    if ~isequal(available_locations, temporal_locations)
        error(['Unexpected temporal locations in %s. Expected [%s], ', ...
            'found [%s].'], input_file, ...
            num2str(temporal_locations), num2str(available_locations));
    end

    exclusions = string(excluded_subjects{group_index});
    if ~isempty(exclusions)
        keep_rows = ~ismember(subject_values, exclusions);
        subject_values = subject_values(keep_rows);
        location_values = location_values(keep_rows);
        error_values = error_values(keep_rows);
    end

    subjects = sort(unique(subject_values));
    fprintf('\n%s group: fitting %d participants\n', ...
        group_label, numel(subjects));

    for subject_index = 1:numel(subjects)
        subject_id = subjects(subject_index);

        for location_index = 1:numel(temporal_locations)
            temporal_location = temporal_locations(location_index);
            row_selector = subject_values == subject_id & ...
                location_values == temporal_location;

            signed_errors = error_values(row_selector);
            signed_errors = signed_errors(isfinite(signed_errors));

            if numel(signed_errors) < minimum_trials_per_fit
                error(['%s (%s) has only %d valid trials at the %d%% ', ...
                    'temporal location; at least %d are required.'], ...
                    subject_id, group_label, numel(signed_errors), ...
                    temporal_location, minimum_trials_per_fit);
            end

            if any(abs(signed_errors) > 80)
                error(['%s (%s) has signed errors outside [-80, 80] at ', ...
                    'the %d%% temporal location.'], subject_id, ...
                    group_label, temporal_location);
            end

            transformed_errors = asind(signed_errors / 80);
            mem_data = struct();
            mem_data.errors = transformed_errors(:)';

            try
                fit_result = MemFit(mem_data, model, 'Verbosity', 0);
                fitted_parameters = fit_result.maxPosterior;
            catch fit_exception
                error(['MemToolbox fitting failed for %s (%s), %d%%: ', ...
                    '%s'], subject_id, group_label, temporal_location, ...
                    fit_exception.message);
            end

            if any(~isfinite(fitted_parameters))
                error(['MemToolbox returned non-finite parameters for ', ...
                    '%s (%s), %d%%.'], subject_id, group_label, ...
                    temporal_location);
            end

            fit_rows(end + 1, :) = {char(subject_id), group_label, ...
                temporal_location, fitted_parameters(mu_index), ...
                fitted_parameters(sigma_index), ...
                fitted_parameters(g_index)}; %#ok<SAGROW>

            fprintf('  %s, %d%%: mu = %.6f, sigma = %.6f, g = %.6f\n', ...
                subject_id, temporal_location, ...
                fitted_parameters(mu_index), ...
                fitted_parameters(sigma_index), ...
                fitted_parameters(g_index));
        end
    end
end

%% Save fitted parameters without modifying the published parameter file
parameter_table = table( ...
    string(fit_rows(:, 1)), ...
    string(fit_rows(:, 2)), ...
    cell2mat(fit_rows(:, 3)), ...
    cell2mat(fit_rows(:, 4)), ...
    cell2mat(fit_rows(:, 5)), ...
    cell2mat(fit_rows(:, 6)), ...
    'VariableNames', ...
    {'subject', 'group', 'temporal_location_pct', 'mu', 'sigma', 'g'} ...
);

output_file = fullfile(output_dir, ...
    'memtoolbox_parameters_reproduced.csv');
writetable(parameter_table, output_file);

fprintf('\nSaved %d fitted parameter rows to:\n%s\n', ...
    height(parameter_table), output_file);
