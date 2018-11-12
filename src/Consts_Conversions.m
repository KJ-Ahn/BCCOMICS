%% Constants in cgs and conversion factors

G        = 6.67259e-8;     %% graviational const.
mH       = 1.6733e-24;     %% gram
kb       = 1.380658e-16;   %% Boltzmann const.
pc       = 3.086e18;       %% pc in cm
Mpc      = pc * 1e6;       %% Mpc in cm
cm_inMpc = 1/Mpc;          %% cm in Mpc
km_inMpc = cm_inMpc * 1e5; %% km in Mpc
yr       = 365.25*24*3600; %% year in sec
Myr      = yr * 1e6;       %% Myr in sec
s_inMyr  = 1/Myr;          %% sec in Myr
c_inkms  = 2.99892458e5;   %% speed of light in km/s
%% multiply to convert Mpc/Myr velocity to km/s velocity
MpcMyr_2_kms = Mpc/Myr /1e5;

