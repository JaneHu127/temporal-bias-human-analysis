%-----------------------------------------------------------------------
% Job saved on 30-Jan-2024 01:35:32 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------


script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));

bids_root = fullfile(repo_root, 'data', 'fMRI', ...
    'subjects_fmridata_BIDS');

output_root = fullfile(repo_root, 'results', 'fMRI', ...
    'first_level_parametric_modulation');
output_root = fullfile(bids_root, 'derivatives', ...
    'spm12', 'first_level', 'temporallocation');

% Subjects with four complete runs
ID = {'010','019'};

spm('defaults', 'fmri');
spm_jobman('initcfg');

for i = 1:numel(ID)
    outputdir = fullfile(output_root, ['sub-' ID{i}]);

    if ~isfolder(outputdir)
        mkdir(outputdir);
    end

    matlabbatch = preproc(ID{i}, bids_root, outputdir);
    spm_jobman('run', matlabbatch);
end
  
  
function matlabbatch = preproc(ID, bids_root, outputdir)
%%
% All runs are stored in the same BIDS functional directory
func_dir = fullfile(bids_root, ['sub-' ID], 'func');

% Select preprocessed 3D images separately for each run
files01 = dir(fullfile(func_dir, ...
    ['swrasub-' ID '_task-temporallocationestimate_run-1*.nii']));
files02 = dir(fullfile(func_dir, ...
    ['swrasub-' ID '_task-temporallocationestimate_run-2*.nii']));
files03 = dir(fullfile(func_dir, ...
    ['swrasub-' ID '_task-temporallocationestimate_run-3*.nii']));


if isempty(files01) || isempty(files02) || isempty(files03)
    error('Missing preprocessed images for sub-%s.', ID);
end

f1 = fullfile({files01.folder}, {files01.name})';
f2 = fullfile({files02.folder}, {files02.name})';
f3 = fullfile({files03.folder}, {files03.name})';


    tsv01 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-1_events.tsv']);
    tsv02 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-2_events.tsv']);
    tsv03 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-3_events.tsv']);
   
    onsets_run1 = tdfread(tsv01);
    onsets_run2 = tdfread(tsv02);
    onsets_run3 = tdfread(tsv03);
  

    rp01 = find_rp_file(func_dir, ID, 1, numel(f1));
    rp02 = find_rp_file(func_dir, ID, 2, numel(f2));
    rp03 = find_rp_file(func_dir, ID, 3, numel(f3));
    


matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 36;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 17; %round(0.880/bin_duration)+1=17,bin_duration=2/36
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = f1;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.name = 'run1';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.onset = onsets_run1.onset;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.tmod = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(1).name = 'error';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(1).param = onsets_run1.ERROR_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(1).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(2).name = 'retention_delay';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(2).param = onsets_run1.ts_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(2).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(3).name = 'sim_triplet_avg';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(3).param = onsets_run1.sim_triplet_avg_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(3).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(4).name = 'signed_DistanceToBoundary_std';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(4).param = onsets_run1.signed_DistanceToBoundary_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod(4).poly = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.orth = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {rp01};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
%%

matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = f2;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.name = 'run2';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.onset = onsets_run2.onset;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.tmod = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(1).name = 'error';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(1).param = onsets_run2.ERROR_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(1).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(2).name = 'retention_delay';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(2).param = onsets_run2.ts_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(2).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(3).name = 'sim_triplet_avg';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(3).param = onsets_run2.sim_triplet_avg_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(3).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(4).name = 'signed_DistanceToBoundary_std';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(4).param = onsets_run2.signed_DistanceToBoundary_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.pmod(4).poly = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond.orth = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg ={rp02};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
%%

matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = f3;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.name = 'run3';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.onset = onsets_run3.onset;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.tmod = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(1).name = 'error';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(1).param = onsets_run3.ERROR_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(1).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(2).name = 'retention_delay';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(2).param = onsets_run3.ts_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(2).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(3).name = 'sim_triplet_avg';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(3).param = onsets_run3.sim_triplet_avg_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(3).poly = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(4).name = 'signed_DistanceToBoundary_std';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(4).param = onsets_run3.signed_DistanceToBoundary_std;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.pmod(4).poly = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond.orth = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {rp03};
matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128;

%%
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'error';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'retention_delay';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 0 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'sim_triplet_avg';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 1 0];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'signed_DistanceToBoundary';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.delete = 0;
end
function rp_file = find_rp_file(func_dir, ID, run_number, number_of_scans)
% Find and validate the motion-parameter file for one run

rp_pattern = sprintf( ...
    'rp_*sub-%s_task-temporallocationestimate_run-%d*.txt', ...
    ID, run_number);

rp_files = dir(fullfile(func_dir, rp_pattern));

if isempty(rp_files)
    error('No motion-parameter file found for sub-%s run-%d.', ...
        ID, run_number);
end

if numel(rp_files) > 1
    error('Multiple motion-parameter files found for sub-%s run-%d.', ...
        ID, run_number);
end

rp_file = fullfile(rp_files(1).folder, rp_files(1).name);
motion_parameters = load(rp_file);

if size(motion_parameters, 2) ~= 6
    error('Motion-parameter file must contain six columns: %s', ...
        rp_file);
end

if size(motion_parameters, 1) ~= number_of_scans
    error(['Motion/scan count mismatch for sub-%s run-%d: ' ...
           '%d motion rows but %d images.'], ...
          ID, run_number, size(motion_parameters, 1), ...
          number_of_scans);
end
end