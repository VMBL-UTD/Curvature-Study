# Arterial Curvature Simulations with MATLAB

## Overview
This repository contains the scripts used to generate idealized arterial geometries with varying plaque size, arterial curvature magnitude, and plaque location relative to the artery's curved centerline. The scripts were developed in the Vascular Mechanobiology Lab (VMBL) at UT Dallas to investigate how arterial geometry affects wall stresses and strains, and subsequently plaque rupture risk.

## Pre-installation
- MATLAB (developed and tested on R2022a)
- [GIBBON toolbox](https://www.gibboncode.org/Installation/) — must be installed before running any scripts in this repository
- FEBio (for finite element solving; see `runSimsHPC.m` for solver calls)

## Usage
The pipeline runs in four sequential steps. Each step depends on the outputs of the previous one, so they should be run in order.

### Step 1: Mesh Generation
In the `main/` folder, run `generateMeshesHPC.m` to generate the 96 meshes, each corresponding to a unique combination of arterial geometry parameters (stenosis severity, curvature magnitude, and curvature direction). This script is designed to run on the UT Dallas HPC cluster but can be adapted for local execution by modifying the relevant job-submission settings.

### Step 2: Running Simulations
In the `main/` folder, run `runSimsHPC.m` to run the FEA solver on each mesh generated in Step 1. As with Step 1, this script targets the UT Dallas HPC cluster but can be tweaked to run locally.

### Step 3: Post-Processing
In the `post processing/` folder, run `postProcessing.m` to compute the mean and max, tensile and compressive, stress and strain values for each simulation, and to format the results for plotting.

### Step 4: Plotting
In the `post processing/` folder:
- Run `boxPlots.m` to generate summary box plots for each measurement.
- Run `spearmanHeatmaps.m` to generate heatmaps visualizing the Spearman rank correlation between arterial geometry parameters and stress/strain outputs.

## Data Usage
Example `.csv` datasets, the Cauchy stress and Green-Lagrange strain tensor components exported from FEBio simulations, are included in the Supplementary Information of the associated manuscript.

[Manuscript citation]
