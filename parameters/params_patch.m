%% Initial condition configuration inside a chosen patch
Lbox_p    = 1/h;         %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p   = 512;         %% # of cells & particles along one axis: make it an even number

ICdir     = '../ICs';  %% directory to place initial condition outputs

oldseedflag    = true;   %% true 
recordseedflag = true;   %% record newly generated seed