%-----------------------------------------------------------------------
% SPM12 Batch Preprocessing - BIDS input
%-----------------------------------------------------------------------

% Resolve the dataset path relative to this script.
% Expected layout:
%   repository_root/
%       spm12_preprocessing_bids.m
%       subjects_fmridata_BIDS/
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
bids_root = fullfile(script_dir, 'subjects_fmridata_BIDS');

% Subject IDs without the "sub-" prefix
IDs = {'001','002','004','005','008','009','010','011','012', ...
       '013','014','015','017','018','019','020','021','022'};

% Initialize SPM
spm('defaults', 'fmri');
spm_jobman('initcfg');

for i = 1:numel(IDs)
    subj_id = IDs{i};
    subj_label = ['sub-' subj_id];
    subj_dir = fullfile(bids_root, subj_label);
    func_dir = fullfile(subj_dir, 'func');
    anat_dir = fullfile(subj_dir, 'anat');

    if ~isfolder(func_dir)
        error('The functional directory for %s does not exist: %s', ...
              subj_label, func_dir);
    end

    % --- Detect the available runs from the 3D functional filenames ---
    run_candidates = dir(fullfile(func_dir, sprintf( ...
        '%s_task-temporallocationestimate_run-*.nii', subj_label)));

    run_numbers = [];
    for f = 1:numel(run_candidates)
        token = regexp(run_candidates(f).name, '_run-(\d+).*\.nii$', ...
                       'tokens', 'once');
        if ~isempty(token)
            run_numbers(end + 1) = str2double(token{1}); %#ok<SAGROW>
        end
    end

    run_numbers = unique(run_numbers);
    run_numbers = run_numbers(:)';

    if isempty(run_numbers)
        error('No functional runs were detected for %s in: %s', ...
              subj_label, func_dir);
    end

    expected_run_numbers = 1:numel(run_numbers);
    if ~isequal(run_numbers, expected_run_numbers)
        error(['The run numbers for %s must start at 1 and be consecutive. ' ...
               'Detected runs: %s'], ...
              subj_label, mat2str(run_numbers));
    end

    nRuns = numel(run_numbers);
    fprintf('%s: detected %d functional runs (%s).\n', ...
            subj_label, nRuns, mat2str(run_numbers));
    scans = cell(1, nRuns);

    for r = 1:nRuns
        % This pattern accepts 3D files such as:
        % sub-001_task-temporallocationestimate_run-1_bold_0001.nii
        % It deliberately excludes SPM outputs whose names start with a prefix.
        file_pattern = sprintf('^%s_task-temporallocationestimate_run-%d.*\.nii$', ...
                               subj_label, r);
        files = spm_select('FPList', func_dir, file_pattern);

        if isempty(files)
            error('No 3D NIfTI files were found for %s run-%d in: %s', ...
                  subj_label, r, func_dir);
        end

        scans{r} = cellstr(strcat(files, ',1'));
    end

    % --- Get the anatomical image ---
    % Each anat directory must contain exactly one uncompressed *_T1w.nii.
    if ~isfolder(anat_dir)
        error('The anatomical directory for %s does not exist: %s', ...
              subj_label, anat_dir);
    end

    anat_nii_list = dir(fullfile(anat_dir, '*_T1w.nii'));
    if numel(anat_nii_list) ~= 1
        error('%s must have exactly one *_T1w.nii file; found %d in: %s', ...
              subj_label, numel(anat_nii_list), anat_dir);
    end

    anat_path = fullfile(anat_nii_list(1).folder, anat_nii_list(1).name);
    anat_file = {[anat_path, ',1']};

    % Use the active SPM installation instead of a hard-coded toolbox path
    tpm_path = fullfile(spm('Dir'), 'tpm', 'TPM.nii');
    if ~isfile(tpm_path)
        error('The SPM TPM file does not exist: %s', tpm_path);
    end

    % -------------------------------------------------------------------
    % Build the preprocessing batch
    % -------------------------------------------------------------------
    matlabbatch = {};

    %% 1. Slice Timing
    matlabbatch{1}.spm.temporal.st.scans = scans;
    matlabbatch{1}.spm.temporal.st.nslices = 36;
    matlabbatch{1}.spm.temporal.st.tr = 2;
    matlabbatch{1}.spm.temporal.st.ta = 0;
    % dcm2niix stores SliceTiming in seconds. This SPM12 configuration
    % expects acquisition times in milliseconds, so multiply by 1000.
    matlabbatch{1}.spm.temporal.st.so = ...
        1000 * [0.88 0 1.1 0.22 1.32 0.44 1.54 0.66 1.76 ...
                0.88 0 1.1 0.22 1.32 0.44 1.54 0.66 1.76 ...
                0.88 0 1.1 0.22 1.32 0.44 1.54 0.66 1.76 ...
                0.88 0 1.1 0.22 1.32 0.44 1.54 0.66 1.76];
    matlabbatch{1}.spm.temporal.st.refslice = 880;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';

    %% 2. Realign
    for r = 1:nRuns
        matlabbatch{2}.spm.spatial.realign.estwrite.data{r}(1) = ...
            cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)', r), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{r}, '.','files'));
    end
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

    %% 3. Coregister
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = anat_file;
    matlabbatch{3}.spm.spatial.coreg.estimate.source(1) = ...
        cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, ...
                  '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','rmean'));
    for r = 1:nRuns
        matlabbatch{3}.spm.spatial.coreg.estimate.other(r) = ...
            cfg_dep(sprintf('Realign: Estimate & Reslice: Resliced Images (Sess %d)', r), ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, ...
                      '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sess', '()',{r}, '.','rfiles'));
    end
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 ...
         0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    %% 4. Normalise
    matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = anat_file;
    for r = 1:nRuns
        matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample(r) = ...
            cfg_dep(sprintf('Realign: Estimate & Reslice: Resliced Images (Sess %d)', r), ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, ...
                      '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sess', '()',{r}, '.','rfiles'));
    end
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpm_path};
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.bb = ...
        [-78 -112 -70; 78 76 85];
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.vox = [3 3 3];
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

    %% 5. Smooth
    matlabbatch{5}.spm.spatial.smooth.data(1) = ...
        cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', ...
        substruct('.','val', '{}',{4}, '.','val', '{}',{1}, ...
                  '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('()',{1}, '.','files'));
    matlabbatch{5}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{5}.spm.spatial.smooth.dtype = 0;
    matlabbatch{5}.spm.spatial.smooth.im = 0;
    matlabbatch{5}.spm.spatial.smooth.prefix = 's';

    % Run the batch for the current subject
    fprintf('\n========== Processing %s ==========\n', subj_label);
    spm_jobman('run', matlabbatch);
end

