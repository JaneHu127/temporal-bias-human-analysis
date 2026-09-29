%% Second-level one-sample t-test: (20% + 40%) > (60% + 80%)
%
% First-level input:
%   con_0005.nii, created by the T contrast [1 1 -1 -1]
%   in location1stlevel.m.
%
% The second-level design has one regressor (the group mean), so its
% positive T contrast is [1]. Do not enter spmT_0005.nii here.

clear;
clc;

%% Configuration
% Resolve all paths relative to this script
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

bids_root = fullfile(script_dir, 'subjects_fmridata_BIDS');

first_level_root = fullfile(bids_root, 'derivatives', ...
    'spm12', 'first_level', 'temporallocation');

second_level_dir = fullfile(bids_root, 'derivatives', ...
    'spm12', 'second_level', 'temporallocation', ...
    '20plus40_gt_60plus80');

subject_ids = {'001','002','004','005','008','009','010','011','012', ...
               '013','014','015','017','018','019','020','021','022'};

first_level_contrast_number = 5;
first_level_contrast_file = sprintf('con_%04d.nii', ...
    first_level_contrast_number);

if exist('spm', 'file') ~= 2 || exist('spm_jobman', 'file') ~= 2
    error('SPM12 is not on the MATLAB path.');
end

if ~isfolder(first_level_root)
    error('First-level directory not found: %s', first_level_root);
end

if ~isfolder(second_level_dir)
    [ok, message] = mkdir(second_level_dir);
    if ~ok
        error('Could not create %s: %s', second_level_dir, message);
    end
end

second_level_spm = fullfile(second_level_dir, 'SPM.mat');
if isfile(second_level_spm)
    error('Second-level SPM.mat already exists: %s', second_level_spm);
end

%% Collect one first-level contrast image per participant
scans = cell(numel(subject_ids), 1);
for i = 1:numel(subject_ids)
    contrast_path = fullfile(first_level_root, ...
        ['sub-' subject_ids{i}], first_level_contrast_file);

    if ~isfile(contrast_path)
        error('Missing first-level contrast for sub-%s: %s', ...
            subject_ids{i}, contrast_path);
    end

    scans{i} = [contrast_path ',1'];
end

%% Specify, estimate, and contrast the second-level model
spm('defaults', 'fmri');
spm_jobman('initcfg');

matlabbatch = {};

% One-sample t-test across participants
matlabbatch{1}.spm.stats.factorial_design.dir = {second_level_dir};
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

% Estimate model
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep( ...
    'Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
              '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Positive group effect: (20% + 40%) > (60% + 80%)
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep( ...
    'Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, ...
              '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ...
    '(20%+40%) > (60%+80%)';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

save(fullfile(second_level_dir, 'second_level_batch.mat'), ...
    'matlabbatch', 'subject_ids', 'scans', ...
    'first_level_contrast_number');

spm_jobman('run', matlabbatch);

fprintf('\nFinished second-level analysis.\n');
fprintf('First-level input: con_%04d.nii ([1 1 -1 -1])\n', ...
    first_level_contrast_number);
fprintf('Output directory: %s\n', second_level_dir);

