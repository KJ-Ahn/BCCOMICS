%% Sanity check for power spectrum normalization
%% Make As, ns, k0, h, omhb2, omch2 consistent with the input (e.g. params.ini) 

scalar_amp(1)             = 2.2e-9
scalar_spectral_index(1)  = 0.961
pivot_scalar              = 0.05
hubble                    = 70.3
ombh2                     = 0.022239
omch2                     = 0.114160

As = scalar_amp 
ns = scalar_spectral_index 
k0 = pivot_scalar
h  = hubble/100
fb = ombh2/(ombh2+omch2)
fc = omch2/(ombh2+omch2)

Pk_zi = load('bccomics_matterpower_z1000.dat');
kh_zi = Pk_zi(:,1); %% in h Mpc^-1 unit
k_zi  = kh_zi*h;    %% in Mpc^-1 unit
Pkm_zi = Pk_zi(:,2)*h^-3; %% into Mpc^3 unit

TF_zi = load('bccomics_transfer_z1000_out.dat');
kh    = TF_zi(:,1); %% in h Mpc^-1 unit
k     = kh*h;       %% in Mpc^-1 unit
col2  = TF_zi(:,2);
col3  = TF_zi(:,3);

TFm   = fc*col2+fb*col3;

Pkm   = As*(k/k0).^(ns-1).*k*2*pi^2.*TFm.^2;
loglog(k, Pkm, k_zi, Pkm_zi) %% plot k(Mpc^-1) vs. Pk(Mpc^3)

