# Optimal Inflation Target with Expectations Driven Liquidity Traps

## System Requirements
All codes are run in MATLAB and heavily utilize fsovle in the optimization toolbox. Pleae make sure this is installed.

Additionally, for the AR(1) portion of the code, mex functions are utilized for speed pickups. The mex functions call interpolation codes written in Fortran. To sucessfully run, make sure you have installed a compatable MATLAB mex Fortran compiler (see MATLAB documentation for a list) and have taken the appropriate steps in MATLAB to run mex functions.

## Folders
There are three folders that contain relevant information
- **AnalyticalJacobianDerivation/**
- **Data/**
- **Figs/**

### AnalyticalJacobianDerivation/
Documentation on the derivation of the analytical Jacobian used in the AR(1) portion of the model

### Data/
Contains excel files and documentation to data sources on Japanese data related to inflation, output gap, and interest rates

### Figs/
Contains all codes used in the paper, organized by the figure they generate.
- **Calib/** contains codes used to calibrate the model used in the paper. Description on calibration method can be found in the paper.
- **Fig_AR1/PiTarg_pT** contains codes used to generate Figure 9 of current paper.
- **Fig_AR1/PFs_and_IRFs** contains codes used to generate Figure 18 of current paper.
- **Fig_AS_AD/** contains codes to generate Figures 3, 12, 13, 14, 16, and 17 of current paper.
- **Fig_EqExistance/** contains codes used to generate figure 10 of current paper.
- **Fig_FR/** contains codes used to generate figure 1 of current paper.
- **Fig_IRFDemand/** contains codes used to generate figure 2 of current paper.
- **Fig_IRFSun/** contains codes used to generate figure 5 of current paper.
- **Fig_LogLin/** contains codes used to generate figure 15 of current paper.
- **Fig_PiTargT/** contains codes used to generate figure 6 of current paper.
- **Fig_PiTargpSA/** contains codes used to generate figure 8 of current paper.
- **Fig_UncProb/** contains codes used to generate figure 7 of current paper.
- **Fig_Welfare/** contains codes used to generate figure 4 of current paper.
- **Final/** a folder that contains all the 'final' versions of figures used in the paper. The .tex file pulls fig files form this folder.
