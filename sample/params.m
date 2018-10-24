%% parameter file for bccomics_setup.m

%%%% Example definitions if the package directory is ~/BCCOMICS, and 
%%%% input files are under ~/BCCOMICS/sample.
%%pkgdir     = '~/BCCOMICS';  %% BCCOMICS directory
%%setupdir   = '~/BCCOMICS/setup_output';  %% directory to record output of bccomics_setup
%%zxestr     = '~/BCCOMICS/sample/output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
%%Cosmology  = '~/BCCOMICS/sample/LCDM.m'; %% '*.m' file containing cosmological parameters that were used as inputs for CAMB
%%TFdir      = '~/BCCOMICS/sample/CAMB_for_mode_finding'; %% directory of CAMB TF outputs

pkgdir     = '..';  %% BCCOMICS directory
setupdir   = '../setup_output';  %% directory to record output of bccomics_setup
zxestr     = 'output_recfast'; %% redshift-(ionized fraction) data; recfast preferred
Cosmology  = 'LCDM.m'; %% '*.m' file containing cosmological parameters that were used as inputs for CAMB
TFdir      = 'CAMB_for_mode_finding'; %% directory of CAMB TF outputs
TFstr1     = [TFdir '/' 'bccomics_transfer_z']; %% CAMB TF output string (head)
TFstr2     = '_out.dat';                        %% CAMB TF output string (tail)
TFzred     = [TFdir '/' 'redshifts.dat'];       %% CAMB TF output redshift file: at least 1000 and 800 should exist.


plotflag   = false     %% true if plots wanted, or false
THflag     = false     %% true if formalism by Tseliakhovich & Hirata and its output are ALSO wanted
OWRTflag   = true      %% true if to overwrite existing data output (under outputdir/deltas and outputdir/deltasTH) wanted.

zzend      = 200;    %% Redshift at which you want to have initial condition.

%% Box configuration for patches. Ncell being an odd number makes FFTing
%% intuitively easier because k runs from -floor(Ncell/2) to floor(Ncell/2),
%% in a symmetric way.
Lbox    = 604;         %% in Mpc unit; let it be (odd number)*4
Ncell   = 151;         %% # of cells along one axis: make it an odd number

%% If you know which cell to choose already, turn on the following flag
%% Useful if bccomics_setup has not finished calculation and you want
%% to rerun bccomics_setup
patchidxinput_flag = false;  %% if true, idx1 & idx2 & idx3 need to be specified
idx1 = 110; %% x-index of your patch (ignored when patchidxinput_flag = false)
idx2 = 24;  %% y-index of your patch (ignored when patchidxinput_flag = false)
idx3 = 115; %% z-index of your patch (ignored when patchidxinput_flag = false)
