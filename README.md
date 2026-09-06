# MEFMA: A Multi-Scale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration

[![Language](https://img.shields.io/badge/Language-MATLAB-blue.svg)](https://www.mathworks.com/)
[![Status](https://img.shields.io/badge/Status-Published-success.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

This repository contains the core algorithmic source code, simulation data generation scripts, quantitative evaluation scripts, and raw datasets for the paper **"MEFMA: A Multiscale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration"**, published in the *IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing (JSTARS)*.

---

## 📖 About the Project

In conditions of extremely low signal-to-noise ratio and complex decorrelation, the interferometric synthetic aperture radar (InSAR) phase is often severely corrupted by dense noise. Existing filtering methods struggle to achieve an effective balance among noise suppression, high-frequency topological detail preservation, and the avoidance of nonphysical spurious fringes, thereby constraining the accuracy of subsequent phase unwrapping.

To address this issue, this article proposes a multiscale empirical fringe-manifold adaptive filtering (MEFMA) framework. The algorithm constructs a V-cycle-like phase decomposition and cascaded reconstruction framework, which is primarily driven alternately by two nonlinear operators:
* **Sparsity-Driven Frequency-Adaptive Operator:** Enables the data-driven extraction of local fringe features based on Hoyer's sparsity, bypassing rigid hyperparameter constraints.
* **Complex Gradient Constant Phase Reconstruction (CGCPR):** A terminal smoothing operator that effectively filters out broadband noise while preserving high-frequency deformation gradients without relying on heavy subspace eigen-decomposition.

Furthermore, the algorithm introduces a coherence-based cascaded weighting and terminal gradient smoothing mechanism to suppress the nonlinear propagation and accumulation of errors.

---

## 📂 Repository Contents

The current repository provides all necessary scripts and data to reproduce the simulation scenarios and evaluate the phase filtering metrics discussed in the manuscript.

### 1. Core Algorithm
* **`MEFMA_Reconstruction.m`**: The primary algorithmic implementation of the MEFMA framework. It executes a rigorous, parameter-free V-cycle multiscale decomposition and cascaded reconstruction. The function is structured into four key processing stages:
  1. **Macroscale Multi-Frequency Phase Decomposition:** Applies adaptive Goldstein frequency stripping across varying window scales.
  2. **Micro-Feature Residual Denoising:** Employs scale-bound fringe manifold projection and coherence-guided phase attenuation to prevent nonlinear noise amplification.
  3. **Coherent Cascaded Reconstruction:** Synthesizes the global phase from coarse to fine scales.
  4. **Empirical Fringe-Manifold Reverse Polishing:** Implements block-based adaptive eigen-fringe Wiener shrinkage, driven by Wax-Kailath Minimum Description Length (MDL) criteria for intrinsic dimensionality estimation. This step ensures extreme preservation of fringe skeletons and steep deformation gradients while eliminating granular noise.

### 2. Simulated Data Generation
We provide the exact MATLAB scripts used to mathematically model the decorrelation noise and generate the complex interferometric fringes.
* `simulate_Cone_Ramp_Peaks.m`: Generates the continuous low-SNR gradient datasets, including the isotropic gradual deformation (Cone), multi-oriented sector radial field (Ramp), and multi-peak topographic phase (Peaks).
* `simulate_Fractal_Composite.m`: Generates the highly non-stationary composite deformation field superimposed on a fractal base.

### 3. Quantitative Evaluation Metrics
These scripts are used to calculate the 6 objective metrics presented in the paper. They are provided to ensure strict experimental reproducibility without ambiguity caused by different programmatic implementations.
* `evaluate_SSIM_RMSE_PSNR.m`: Evaluates Structural Similarity (SSIM), Root Mean Square Error (RMSE), and Peak Signal-to-Noise Ratio (PSNR) against the ground-truth phase.
* `evaluate_PSD_SPD_RPN.m`: Evaluates local phase smoothness via Phase Standard Deviation (PSD) and Sum of Absolute Phase Differences (SPD), and calculates the Residue Point Number (RPN) for topological continuity.

### 4. Datasets
* `Simulated_Datasets.mat`: Contains the pre-generated simulated complex phase fields (wrapped noise phase, true phase, and coherence maps).
* `Real_Datasets.mat`: Contains the highly decorrelated real-world spaceborne SAR datasets used in Section IV of the manuscript, including the **Myanmar earthquake**, **Turkey earthquake**, and **Etna volcano**.

---

## 🚀 Usage 

All scripts are written in **MATLAB**. To run the scripts:

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/YourUsername/MEFMA.git
   cd MEFMA
   ```
2. Ensure you have the **Signal Processing Toolbox** and **Image Processing Toolbox** installed in your MATLAB environment.
3. Run `MEFMA_Reconstruction.m` to apply the filtering algorithm to your own wrapped phase data. 
4. Run the evaluation scripts directly to reproduce the paper's quantitative results. The scripts are configured to automatically load the corresponding variables from the `.mat` files provided in this repository.

---

## 📝 Citation

If you find our code, datasets, or evaluation scripts useful for your research, please cite our paper:

**Plain Text:**
> J. Li, H. Fan, H. Zhuang, Y. Yeerdingdala and S. Du, "MEFMA: A Multiscale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 19, pp. 24910-24926, 2026, doi: 10.1109/JSTARS.2026.3715361.

**BibTeX:**
```bibtex
@article{Li2026MEFMA,
  author={Li, Jian and Fan, Hongdong and Zhuang, Huifu and Yeerdingdala, Yeerda and Du, Sen},
  journal={IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing}, 
  title={MEFMA: A Multiscale Empirical Fringe-Manifold Adaptive Filtering Framework for Robust InSAR Phase Restoration}, 
  year={2026},
  volume={19},
  pages={24910-24926},
  doi={10.1109/JSTARS.2026.3715361}
}
```
