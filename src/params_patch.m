%% Initial condition configuration inside a chosen patch
Lbox_p     = 1/h       %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p    = 512       %% # of cells & particles along one axis: make it an even number

plotflag_p = false     %% true if plots of initial conditions wanted

ICdir      = '../ICs';  %% (mother) directory for initial condition output
