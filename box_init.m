%% Box configuration for patches. Ncell being an odd number makes FFTing
%% intuitively easier because k runs from -floor(Ncell/2) to floor(Ncell/2),
%% in a symmetric way.
Vbox    = Lbox^3;      %% box volume in Mpc^3 unit
Lcell   = Lbox/Ncell;  %% cell size in Mpc unit; 4 Mpc is fiducial, for correct DeltaT fit
Vcell   = Lcell^3;     %% cell volume in Mpc^3 unit
Nmode   = Ncell;       %% number of k modes along one axis
Nhalf   = floor(Ncell/2); 
Nc      = Nhalf+1;     %% index for center of k-space 
kunit   = 2*pi/Lbox;   %% unit k in Mpc^-1
