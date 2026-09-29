%-----------------------------------------------------------------------
% Job saved on 30-Jan-2024 01:35:32 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% Resolve all paths relative to this script
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

bids_root = fullfile(script_dir, 'subjects_fmridata_BIDS');
output_root = fullfile(bids_root, 'derivatives', ...
    'spm12', 'first_level', 'temporallocation');

% Subjects with four complete runs
ID = {'001','002','004','005','008','009','011','012', ...
      '013','014','015','017','018','020','021','022'};

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
files04 = dir(fullfile(func_dir, ...
    ['swrasub-' ID '_task-temporallocationestimate_run-4*.nii']));

if isempty(files01) || isempty(files02) || ...
        isempty(files03) || isempty(files04)
    error('Missing preprocessed images for sub-%s.', ID);
end

f1 = fullfile({files01.folder}, {files01.name})';
f2 = fullfile({files02.folder}, {files02.name})';
f3 = fullfile({files03.folder}, {files03.name})';
f4 = fullfile({files04.folder}, {files04.name})';


    tsv01 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-1_events.tsv']);
    tsv02 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-2_events.tsv']);
    tsv03 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-3_events.tsv']);
    tsv04 = fullfile(func_dir, ...
        ['sub-' ID '_task-temporallocationestimate_run-4_events.tsv']);

    rp01 = find_rp_file(func_dir, ID, 1, numel(f1));
    rp02 = find_rp_file(func_dir, ID, 2, numel(f2));
    rp03 = find_rp_file(func_dir, ID, 3, numel(f3));
    rp04 = find_rp_file(func_dir, ID, 4, numel(f4));
    
onsets_run1 = tdfread(tsv01);
matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 36;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 18;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = [f1];
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = '20%';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = onsets_run1.onset(onsets_run1.temporallocation==1);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = '40%';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = onsets_run1.onset(onsets_run1.temporallocation==5)';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).name = '60%';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).onset = onsets_run1.onset(onsets_run1.temporallocation==9);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).name = '80%';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).onset = onsets_run1.onset(onsets_run1.temporallocation==13);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(4).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {rp01};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
%%
onsets_run2 = tdfread(tsv02);
matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = [f2];
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).name = '20%';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).onset = onsets_run2.onset(onsets_run2.temporallocation==2);
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).name = '40%';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).onset = onsets_run2.onset(onsets_run2.temporallocation==6);
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).name = '60%';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).onset = onsets_run2.onset(onsets_run2.temporallocation==10);
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).name = '80%';
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).onset = onsets_run2.onset(onsets_run2.temporallocation==14);
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(4).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg ={rp02};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
%%
onsets_run3 = tdfread(tsv03);
matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = [f3];
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).name = '20%';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).onset = onsets_run3.onset(onsets_run3.temporallocation==3);
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).name = '40%';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).onset = onsets_run3.onset(onsets_run3.temporallocation==7);
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).name = '60%';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).onset = onsets_run3.onset(onsets_run3.temporallocation==11);
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(3).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).name = '80%';
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).onset = onsets_run3.onset(onsets_run3.temporallocation==15);
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(4).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {rp03};
matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128;
%%
onsets_run4 = tdfread(tsv04);
matlabbatch{1}.spm.stats.fmri_spec.sess(4).scans = [f4];
%%
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).name = '20%';
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).onset = onsets_run4.onset(onsets_run4.temporallocation==4);
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).name = '40%';
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).onset = onsets_run4.onset(onsets_run4.temporallocation==8);
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).name = '60%';
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).onset = onsets_run4.onset(onsets_run4.temporallocation==12);
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(3).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).name = '80%';
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).onset = onsets_run4.onset(onsets_run4.temporallocation==16);
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).duration = 10;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond(4).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {rp04};
matlabbatch{1}.spm.stats.fmri_spec.sess(4).hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat ile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '20%';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '40%';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = '60%';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = '80%';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = '1 1 -1 -1';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.delete = 0;
end

function rp_file = find_rp_file(func_dir, ID, run_number, number_of_scans)
% Find the motion-parameter file for one run

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
