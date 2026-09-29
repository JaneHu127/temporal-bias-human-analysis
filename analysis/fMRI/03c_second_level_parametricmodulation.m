%% SPM12 second-level one-sample t-tests for the four pmod effects
%
% Each first-level contrast is the mean effect across that participant's
% available runs. Four-run participants and the two three-run participants
% (sub-010 and sub-019) can therefore enter the same random-effects model.
%
% First-level contrast mapping (verified from every first-level SPM.mat):
%   con_0001.nii = error
%   con_0002.nii = retention delay (ts)
%   con_0003.nii = image similarity (sim_triplet_avg)
%   con_0004.nii = distance to boundary (signed_distance_to_boundary)
%
% A separate one-sample t-test is estimated for each pmod. Each model gets:
%   Positive group effect: [1]
%   Negative group effect: [-1]
%
% The script deliberately stops if an output SPM.mat already exists so that
% a previous second-level analysis is not silently overwritten.

clear;
clc;

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

%% Paths
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

repo_root = fileparts(fileparts(script_dir));

first_level_root = fullfile(repo_root, 'results', 'fMRI', ...
    'first_level_parametric_modulation');

second_level_root = fullfile(repo_root, 'results', 'fMRI', ...
    'second_level_parametric_modulation');

%% Participants
% This uses all 18 completed first-level analyses. To reproduce the older
% 16-participant sample, remove '010' and '019' from this list.
subject_ids = {'001','002','004','005','008','009','010','011','012', ...
               '013','014','015','017','018','019','020','021','022'};

%% Models
models = struct( ...
    'folder',    {'PM01_error', ...
                  'PM02_retention_delay', ...
                  'PM03_image_similarity', ...
                  'PM04_signed_distance_to_boundary'}, ...
    'label',     {'error', ...
                  'retention_delay', ...
                  'image_similarity', ...
                  'signed_distance_to_boundary'}, ...
    'con_index', {1, 2, 3, 4});

if ~isfolder(first_level_root)
    error('First-level root does not exist: %s', first_level_root);
end

if ~isfolder(second_level_root)
    mkdir(second_level_root);
end

%% Specify, estimate, and contrast each second-level model
for model_index = 1:numel(models)
    model = models(model_index);
    output_dir = fullfile(second_level_root, model.folder);

    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    output_spm = fullfile(output_dir, 'SPM.mat');
    if isfile(output_spm)
        error(['Second-level SPM.mat already exists. Refusing to overwrite: ' ...
               '%s'], output_spm);
    end

    scans = cell(numel(subject_ids), 1);
    for subject_index = 1:numel(subject_ids)
        subject_id = subject_ids{subject_index};
        contrast_file = fullfile( ...
            first_level_root, ...
            ['sub-' subject_id], ...
            sprintf('con_%04d.nii', model.con_index));

        if ~isfile(contrast_file)
            error('Missing first-level contrast: %s', contrast_file);
        end

        scans{subject_index} = [contrast_file ',1'];
    end

    matlabbatch = {};

    % 1. One-sample t-test design
    matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
    matlabbatch{1}.spm.stats.factorial_design.cov = ...
        struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = ...
        struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    % 2. Classical model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep( ...
        'Factorial design specification: SPM.mat File', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
                  '.','val', '{}',{1}), ...
        substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % 3. Group-level positive and negative effects
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep( ...
        'Model estimation: SPM.mat File', ...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, ...
                  '.','val', '{}',{1}), ...
        substruct('.','spmmat'));

    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ...
        sprintf('%s_positive', model.label);
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ...
        sprintf('%s_negative', model.label);
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.delete = 0;

    fprintf('\nRunning second-level model %d/%d: %s\n', ...
        model_index, numel(models), model.label);
    spm_jobman('run', matlabbatch);
end

fprintf('\nAll four second-level one-sample models completed.\n');

