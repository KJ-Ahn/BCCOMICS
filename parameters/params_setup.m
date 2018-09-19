%% parameter file for bccomics_setup.m

setupdir   = '../setup_output';  %% directory for output from bccomics_setup
TFdir      = '../sample/CAMB_for_mode_finding'; %% directory of CAMB TF outputs
TFstr1     = [TFdir '/' 'bccomics_transfer_z']; %% CAMB TF output string (head)
TFstr2     = '_out.dat';                        %% CAMB TF output string (tail)
TFzred     = [TFdir '/' 'redshifts.dat'];       %% CAMB TF output redshift file: at least 1000 and 800 should exist.

zxestr     = '../sample/output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
Cosmology  = '../sample/LCDM.m'; %% '*.m' file containing cosmological parameters that were used as inputs for CAMB
plotflag   = false     %% true if plots wanted, or false
THflag     = false     %% true if formalism by Tseliakhovich & Hirata and its output are ALSO wanted
OWRTflag   = true      %% true if to overwrite existing data output (under outputdir/deltas and outputdir/deltasTH) wanted.

zzend      = 200;    %% Redshift at which you want to have initial condition.

%% Box configuration for patches. Ncell being an odd number makes FFTing
%% intuitively easier because k runs from -floor(Ncell/2) to floor(Ncell/2),
%% in a symmetric way.
Lbox    = 604;         %% in Mpc unit; let it be (odd number)*4
Ncell   = 151;         %% # of cells along one axis: make it an odd number

%% Initial condition configuration inside a chosen patch
Lbox_p    = 1/h;         %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p   = 512;         %% # of cells & particles along one axis: make it an even number