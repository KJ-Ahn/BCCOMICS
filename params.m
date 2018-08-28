%% parameter file for bccomics
%% DO NOT CHANGE 

outputdir='setup_output';  %% output directory name
TFstr1   ='CAMB_for_mode_finding/bccomics_transfer_z'; %% CAMB TF output string (head)
TFstr2   ='_out.dat';                                  %% CAMB TF output string (tail)
zxestr   ='output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
Cosmology='LCDM.m'; %% '*.m' file containing cosmological parameters
plotflag =false     %% true if plots wanted, or false
matlabflag=false     %% true if using MATLAB, false if using gnu octave
zzend     = 200;    %% Redshift at which you want to have initial condition.

%% Box configuration for patches. Ncell being an odd number makes FFTing
%% intuitively easier because k runs from -floor(Ncell/2) to floor(Ncell/2),
%% in a symmetric way.
Lbox    = 604;         %% in Mpc unit; let it be (odd number)*4
Ncell   = 151;         %% # of cells along one axis: make it an odd number
