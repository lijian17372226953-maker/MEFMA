# MEFMA: A Multi-Scale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration

This repository contains the simulation data generation scripts, quantitative evaluation scripts, and raw datasets for the paper "**MEFMA: A Multi-Scale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration**", currently under review in the *IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing* (JSTARS).

> **⚠️ Important Notice Regarding the Core Algorithm:** 
> To comply with journal policies during the double-blind/peer-review process, the core MEFMA algorithmic source code (`MEFMA.m` and related functional dependencies) is temporarily withheld. **The complete core algorithm will be fully open-sourced in this repository immediately upon the formal acceptance of the manuscript.**

## 📂 Repository Contents

The current repository provides all necessary scripts and data to reproduce the simulation scenarios and evaluate the phase filtering metrics discussed in the manuscript.

### 1. Simulated Data Generation
We provide the exact MATLAB scripts used to mathematically model the decorrelation noise and generate the complex interferometric fringes.
* `simulate_Cone_Ramp_Peaks.m`: Generates the continuous low-SNR gradient datasets, including the isotropic gradual deformation (Cone), multi-oriented sector radial field (Ramp), and multi-peak topographic phase (Peaks).
* `simulate_Fractal_Composite.m`: Generates the highly non-stationary composite deformation field superimposed on a fractal base.

### 2. Quantitative Evaluation Metrics
These scripts are used to calculate the 6 objective metrics presented in the paper. They are provided to ensure strict experimental reproducibility without ambiguity caused by different programmatic implementations.
* `evaluate_SSIM_RMSE_PSNR.m`: Evaluates Structural Similarity (SSIM), Root Mean Square Error (RMSE), and Peak Signal-to-Noise Ratio (PSNR) against the ground-truth phase.
* `evaluate_PSD_SPD_RPN.m`: Evaluates local phase smoothness via Phase Standard Deviation (PSD) and Sum of Absolute Phase Differences (SPD), and calculates the Residue Point Number (RPN) for topological continuity.

### 3. Datasets
* `Simulated_Datasets.mat`: Contains the pre-generated simulated complex phase fields (wrapped noise phase, true phase, and coherence maps).
* `Real_Datasets.mat`: Contains the highly decorrelated real-world spaceborne SAR datasets used in Section IV of the manuscript, including the **Myanmar earthquake**, **Turkey earthquake**, and **Etna volcano**.

## 🚀 Usage 

All scripts are written in **MATLAB**. To run the scripts:
1. Clone this repository to your local machine.
2. Ensure you have the Signal Processing Toolbox and Image Processing Toolbox installed in your MATLAB environment.
3. Run the evaluation scripts directly. The scripts are configured to automatically load the corresponding variables from the `.mat` files provided in this repository.

