%% parameter file for bccomics_setup.m

%%pkgdir     = '~/Documents/PClinuxOSDATA/BCCOMICS_for_release/BCCOMICS';  %% BCCOMICS directory
%%setupdir   = '~/Documents/PClinuxOSDATA/BCCOMICS_for_release/BCCOMICS/setup_output';  %% directory to record output of bccomics_setup
%%TFdir      = '~/Documents/PClinuxOSDATA/BCCOMICS_for_release/BCCOMICS/sample/CAMB_for_mode_finding'; %% directory of CAMB TF outputs
%%zxestr     = '~/Documents/PClinuxOSDATA/BCCOMICS_for_release/BCCOMICS/sample/output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
%%Cosmology  = '~/Documents/PClinuxOSDATA/BCCOMICS_for_release/BCCOMICS/sample/LCDM.m'; %% '*.m' file containing cosmological parameters that were used as inputs for CAMB

pkgdir     = '..';  %% BCCOMICS directory
setupdir   = '../setup_output';  %% directory to record output of bccomics_setup
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
