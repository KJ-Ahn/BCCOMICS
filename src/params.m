%% parameter file for bccomics.m
%% "Cosmology" and "zzend" should be inherited from params_setup.m.

Cosmology = 'LCDM.m';           %% Cosmology in params_setup.m
setupdir  = '../setup_output';  %% outputdir in params_setup.m
ICdir     = '../ICs';           %% dir to output initial conditions

plotflag  = false     %% true if plots wanted, or false
THflag    = false     %% true if formalism by Tseliakhovich & Hirata and its output are ALSO wanted
OWRTflag  = true      %% true if to overwrite existing data output (under outputdir/deltas and outputdir/deltasTH) wanted.

zzend     = 200;    %% Redshift at which you want to have initial condition.



%% Configuration for patch (as a new box). Ncell_p need to be even numbered.
%% enzo, for example, does not allow odd numbered grid.
%% Lbox_p should be smaller than Lcell in box_init.m.
%% Grid resolution = # of CDM particles = Ncell_p^3

Lbox_p    = 1/h;         %% in Mpc unit; let it be (odd number)*4
Ncell_p   = 512;         %% # of cells along one axis: make it an even number
