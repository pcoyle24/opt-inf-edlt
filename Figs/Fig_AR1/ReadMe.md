## Codes to run

### Pitarg_pT Folder
Codes in this folder should be run **first**. Run run_val_pspd.m. This will run two files:
- (1) script_val_mex_pitarg_pspd_grid.m to generate a number of .mat files stored in \Fig_AR1_savedata101
- (2) val_savedata_mex_pspd.m to call data from Fig_AR1_savedata101 generate Ev_OptInf_sun_5bps.eps

### PFs_and_IRFs Folder
**If you have not run the codes in Pitarg_pT and do not have .mat files saved in the savedata_101 folder, do that first**
- run script_val_mex.m to generate Pfs_sunspot.eps
- run plot_irfs.m to generate IRFs_sunspot.eps 

### Accuracy Folder
Run accuracy_mex.m to generate accuracy.txt, which is reported in the appendix of the paper. 