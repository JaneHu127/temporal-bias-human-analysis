# Temporal Location Bias in Episodic Memory

Data and analysis code associated with the study:

> **The human central gyrus tracks asymmetric temporal location bias in episodic memory representation**  
> Xinyue Hu, Xiuru Cai, Qihong Lv, Lei Wang, Xuanlong Zhu, Jianxuan Yang, Huimin Wang, Chunxia Li, and Sze Chai Kwok

## Overview

This repository accompanies a study of how people reconstruct the temporal location of events in episodic memory. Participants viewed a continuous movie and later estimated where individual frames had appeared within temporal intervals. The analyses examine:

- central-tendency and early-versus-late asymmetries in temporal memory errors;
- bias, precision, and guess-rate parameters estimated with MemToolbox;
- trial-level predictors of signed temporal memory error;
- model-free estimates of temporal-location bias; and
- fMRI effects of temporal location and trial-by-trial memory error.

The repository contains de-identified behavioral data and the MATLAB, R, and Python analysis scripts. Raw and preprocessed MRI data will be available separately on OpenNeuro. The dataset accession number and DOI will be added after publication.
images
## Repository structure

```text
.
|-- README.md
|-- data/
|   |-- eventboundary.xlsx
|   |-- behavioral/
|   |   |-- fMRI_group/
|   |   |   |-- participants.csv
|   |   |   |-- fmri_group_trials.csv
|   |   |   `-- individual_data/
|   |   |-- replication_group/
|   |   |   |-- participants.csv
|   |   |   |-- replication_group_trials.csv
|   |   |   `-- individual_data/
|   |   |-- control_group/
|   |   |   |-- participants.csv
|   |   |   |-- control_group_trials.csv
|   |   |   `-- individual_data/
|   |   |-- memtoolbox_parameters.csv
|   |   `-- participant_rt_summary.csv
|   `-- fMRI/                         # MRI data are not distributed here
`-- analysis/
    |-- behavioral/
    |   |-- 01_fit_memtoolbox_parameters.m
    |   |-- 02_behavioral_statistics.R
    |   |-- 03_trial_level_lmm.R
    |   |-- 04_model_free_bias.py
    |   `-- 05_plot_model_free_bias.py
    `-- fMRI/
        |-- 01_preproc.m
        |-- 02a_first_level_glm_location_for16sub.m
        |-- 02b_first_level_glm_location_for2sub.m
        |-- 02c_second_level_onesample_location.m
        |-- 03a_first_level_glm_parametric_modulation_for16sub.m
        |-- 03b_first_level_glm_parametric_modulation_for2sub.m
        `-- 03c_second_level_parametricmodulation.m
```

Analysis outputs are written to a locally created `results/` directory and are not part of the source data.

## Behavioral datasets

| Sample | Participants | Trial rows | Description |
|---|---:|---:|---|
| fMRI group | 18 | 4,200 | Behavioral data collected during the fMRI experiment. Sixteen participants completed four runs; `sub-10` and `sub-19` completed three runs. |
| Behavioral replication group | 30 | 8,585 | Independent behavioral replication. The `valid_response` column identifies valid and invalid responses. |
| No-encoding control group | 9 | 2,160 | Control experiment performed without prior movie encoding. |

The group-level CSV files contain all trial rows for a sample. The `individual_data/` directories provide the same data separated by participant. Participant tables contain de-identified subject IDs and demographic variables.

### Main behavioral variables

- `subject`, `trial`, `run` or `block`: participant and trial identifiers;
- `interval_width`: length of the tested temporal interval;
- `temporal_location(%)` or `target_location`: target position at 20%, 40%, 60%, or 80% of the interval;
- `reported_location(%)` or `reported_location`: participant response;
- `signed_error`: reported minus target temporal location;
- `absolute_error`: absolute temporal estimation error;
- `response_time`, `confidence`, or related fields: response-time and confidence measures where available;
- `onset`, `retention_delay`, `image_similarity`, and `distance_to_boundary`: trial-level variables used in the fMRI-group analyses.

`memtoolbox_parameters.csv` contains the published participant-by-location estimates of directional bias (`mu`), variability (`sigma`), and guess rate (`g`). The file uses the Greek column names `μ` and `σ`. `eventboundary.xlsx` contains manually identified movie shot-boundary information used to derive distance-to-boundary measures.

## Software requirements

### MATLAB

- MATLAB
- [MemToolbox](https://visionlab.github.io/MemToolbox/) for the behavioral mixture-model analysis
- SPM12 for MRI preprocessing and GLM analyses (the batch files were prepared for SPM12 revision 7771)

### R

The behavioral scripts use the following packages:

```r
readr, dplyr, tidyr, rstatix, afex,
lme4, lmerTest, performance, effectsize, ggplot2
```

### Python

- Python 3
- NumPy
- SciPy
- Matplotlib

For example:

```bash
python -m pip install numpy scipy matplotlib
```

No exact software environment or lock file is currently provided. Package versions should therefore be recorded when reproducing the analyses.

## Behavioral analysis

Run commands from the repository root unless noted otherwise.

### 1. Fit the mixture model

Add MemToolbox to the MATLAB path, then run:

```matlab
run('analysis/behavioral/01_fit_memtoolbox_parameters.m')
```

This fits the bias-aware mixture model separately for each participant and temporal location. It writes:

```text
results/behavioral/memtoolbox_parameters_reproduced.csv
```

The script applies the manuscript-defined exclusion of `sub-12` and `sub-13` from the fMRI sample for the parameter analysis. The published parameter table remains available at `data/behavioral/memtoolbox_parameters.csv`.

### 2. Reproduce the behavioral statistics

```bash
Rscript analysis/behavioral/02_behavioral_statistics.R
```

This script reads the published parameter table and writes statistical summaries and a text report to `results/behavioral/`.

### 3. Fit the trial-level mixed-effects model

```bash
Rscript analysis/behavioral/03_trial_level_lmm.R INPUT.csv OUTPUT_DIRECTORY
```

The input must contain the columns validated near the beginning of the script, including `signed_distance_to_boundary`. That derived column is not present in the distributed `fmri_group_trials.csv`, so it must be prepared before this script is run.

### 4. Compute and plot model-free bias

```bash
python analysis/behavioral/04_model_free_bias.py
python analysis/behavioral/05_plot_model_free_bias.py
```

The first script calculates participant-level and group-level early-versus-late bias measures. The second script reads those outputs and creates PNG and vector PDF figures in `results/behavioral/`.

## fMRI analysis

The MRI scripts are supplied for transparency, but the MRI images and derived SPM files are not distributed in this repository. Before running the scripts, place the local dataset in the expected BIDS-like structure and set `bids_root` consistently in each MATLAB file.

An expected subject layout is:

```text
subjects_fmridata_BIDS/
`-- sub-001/
    |-- anat/
    |   `-- sub-001_T1w.nii
    `-- func/
        |-- sub-001_task-temporallocationestimate_run-1_bold_0001.nii
        |-- ...
        `-- sub-001_task-temporallocationestimate_run-1_events.tsv
```

The preprocessing pipeline in `01_preproc.m` includes slice-timing correction, realignment, coregistration, normalization to MNI space, and 6 mm FWHM spatial smoothing. Acquisition-specific settings in the script correspond to the study sequence: TR = 2 s, 36 slices, and multiband factor 4.

Two first-level GLMs are provided:

1. **Temporal-location GLM (`02a`/`02b`)**: separate 20%, 40%, 60%, and 80% regressors, including the contrast `(20% + 40%) > (60% + 80%)`.
2. **Parametric-modulation GLM (`03a`/`03b`)**: signed memory error, retention delay, image similarity, and distance to boundary as non-orthogonalized parametric modulators.

The `a` scripts analyze the 16 participants with four complete runs. The `b` scripts analyze `sub-010` and `sub-019`, who have three complete runs. The corresponding `c` scripts perform second-level one-sample t-tests.

Some fMRI scripts currently resolve `subjects_fmridata_BIDS` relative to `analysis/fMRI/`, whereas others resolve it under `data/fMRI/`. Because the MRI dataset is not included, users should review and harmonize these path definitions before execution. Run the scripts in numerical order only after confirming the input paths, event files, motion regressors, subject lists, and output directories.

## Related model code

The TCM-based recurrent neural network model is maintained separately:

- <https://github.com/qihongl/temporal-position-bias>

## Citation

If you use these data or scripts, please cite the associated article. Full journal and DOI information will be added when available.

## Data privacy

The behavioral files use coded participant identifiers and do not contain names or direct contact information. Users remain responsible for handling the data in accordance with applicable ethical approvals and data-protection requirements.

## License

No reuse license is currently included. Until a license is added, please contact the authors before redistributing or reusing the data or code beyond citation and inspection.
