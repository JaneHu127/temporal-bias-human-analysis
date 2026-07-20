# Temporal Bias Human Analysis

This repository contains the behavioral data and analysis code used to estimate the MemToolbox bias parameter μ, together with the participant-level MemToolbox results and the R code used to analyze and visualize μ across temporal-location conditions.

## Repository structure

```text
temporal-bias-human-analysis/
├── code/
│   ├── compute_memresult.m
│   └── anovaplot.R
│
├── data/
│   ├── raw_behavioral_data/
│   └── memtoolbox_results/
│
└── README.md

Data
Trial-level behavioral data

The data/signederror_eachsub/ folder contains trial-level signed-error data from 18 participants. 

Signed error was defined as:

signed error (%) = reported temporal position − correct temporal position

These trial-level signed-error data were fitted using MemToolbox to estimate the bias parameter μ separately for each participant and temporal-location condition.

MemToolbox results

The data/memresults.xlsx file contains participant-level parameter estimates obtained from MemToolbox for 16 participants.

For each participant and temporal-location condition, the following parameters were estimated:

μ: systematic directional bias
σ: variability in memory responses
g: probability of random guessing

Participants 12 and 13 were excluded from analyses involving the model-derived parameters because their individual MemToolbox fits failed to converge.


Analysis code
compute_memresult.m

This MATLAB script:

Reads the trial-level signed-error data from 18 participants.
Separates the data by participant and temporal-location condition.
Fits the `Orientation(WithBias(StandardMixtureModel()), [1 3])` model using MemToolbox separately for each participant and temporal-location condition.
Extracts the maximum a posteriori (MAP) estimates of μ, SD, and g from `fitResult.maxPosterior`.
Estimates the parameters μ, σ, and g.
Excludes participants whose individual model fits failed to converge.
Saves the participant-level MemToolbox results to memresults.xlsx.

MemToolbox must be installed and added to the MATLAB path before running this script.

anovaplot.R

This R script:

Reads the participant-level MemToolbox results from memresults.xlsx.
Extracts the μ estimates for the four temporal-location conditions.
Performs the statistical analyses of μ across temporal locations.
Calculates the relevant descriptive statistics.
Visualizes the mean μ estimates across the four temporal-location conditions.

Analysis workflow

Run the scripts in the following order:

Run code/compute_memresult.m in MATLAB.
Confirm that data/memresults.xlsx has been generated successfully.
Run code/anovaplot.R in R to analyze and visualize the μ estimates.

The overall workflow is：

Trial-level signed-error data from 18 participants
                         ↓
             MemToolbox model fitting
                         ↓
            Estimation of μ, σ, and g
                         ↓
 Exclusion of two participants with non-convergent fits
                         ↓
       MemToolbox results from 16 participants
                         ↓
          R analysis and visualization of μ

Software requirements
MATLAB
MATLAB
MemToolbox
R

The required R packages are listed at the beginning of anovaplot.R.

Notes
The temporal-location conditions correspond to 20%, 40%, 60%, and 80% of the temporal interval.
Local file paths may need to be modified before running the scripts on another computer.
