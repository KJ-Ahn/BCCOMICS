%% Cosmology input : Should match CAMB input parameters !!!
ns = 0.961;                      %% scalar_spectral_index(1)
nrun = 0;                        %% scalar_nrun(1)
nrunrun = 0;                     %% scalar_nrunrun(1)
k0 = 0.05;                       %% pivot_scalar
As = 2.2e-9;                     %% scalar_amp(1)
h  = 0.703;                      %% hubble/100
H0 = 100*h*km_inMpc/s_inMyr;     %% H0 in Myr^-1
ombh2 = 0.022239;                %% Omb(0)*h^2
omch2 = 0.114160;                %% Omc(0)*h^2
Om0   = (ombh2+omch2)/h^2;       %% Omega_m(0)
Omgamma    = 4.98569497e-5;      %% Omega_photon(0)
Omneutrino = 3.44215495e-5;      %% Omega_neutrino(0)
Omr0       = Omgamma+Omneutrino; %% Omega_radiation(0) (assuming zero mass for neutrinos)
OmLambda0  = 1-Om0-Omr0;         %% Omega_Lambda(0)
TCMB0 = 2.726;                   %% CMB_temperature(0)
tgamma = 1/8.55e-13 /1e6;        %% in Myr unit (Naoz & Barkana 2005, MNRAS 362, 1047)
