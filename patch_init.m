%% Patch configuration.
Vbox_p  = Lbox_p^3;      %% patch volume in Mpc^3 unit
Lcell_p = Lbox_p/Ncell_p;  %% cell size in Mpc unit
Vcell_p = Lcell_p^3;     %% cell volume in Mpc^3 unit
Nmode_p = Ncell_p;       %% number of k modes along one axis
Nhalf_p = Ncell_p/2; 
Nc_p    = Nhalf_p+1;     %% index for center of k-space 
kunit_p = 2*pi/Lbox_p;   %% unit k in Mpc^-1
