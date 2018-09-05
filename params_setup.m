%% parameter file for bccomics_setup.m

outputdir  = 'setup_output';  %% output directory name
TFstr1     = 'CAMB_for_mode_finding/bccomics_transfer_z'; %% CAMB TF output string (head)
TFstr2     = '_out.dat';                                  %% CAMB TF output string (tail)
zxestr     = 'output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
Cosmology  = 'LCDM.m'; %% '*.m' file containing cosmological parameters that were used as inputs for CAMB
plotflag   = false     %% true if plots wanted, or false
THflag     = false     %% true if formalism by Tseliakhovich & Hirata and its output are ALSO wanted
OWRTflag   = true      %% true if to overwrite existing data output (under outputdir/deltas and outputdir/deltasTH) wanted.

zzend      = 200;    %% Redshift at which you want to have initial condition.

%% Box configuration for patches. Ncell being an odd number makes FFTing
%% intuitively easier because k runs from -floor(Ncell/2) to floor(Ncell/2),
%% in a symmetric way.
Lbox    = 604;         %% in Mpc unit; let it be (odd number)*4
Ncell   = 151;         %% # of cells along one axis: make it an odd number
