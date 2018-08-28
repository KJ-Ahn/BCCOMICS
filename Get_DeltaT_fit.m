%% Prepares for DeltaT fit.
%% Also generates DeltaT at azend.

global beta gamma;
global signDT alpha coeff_Delta_T;

%% construct 3D field of DTA
%% 0.279 is the coefficient, but depending on cosmology this
%% can be somewhat different
aA    = 0.01;
DpgA  = interp1(azz, Dplus_grow   , aA, 'spline');
DpdA  = interp1(azz, Dplus_decay  , aA, 'spline');
DmsA  = interp1(azz, Dminus_stream, aA, 'spline');
DTA3D = 0.279*(Deltagro*DpgA + Deltadec*DpdA -fc*(Deltacom + Deltastr*DmsA));

%% construct 3D field of DTB
%% 0.599 is the coefficient, but depending on cosmology this
%% can be somewhat different
aB    = 0.1;
DpgB  = interp1(azz, Dplus_grow   , aB, 'spline');
DpdB  = interp1(azz, Dplus_decay  , aB, 'spline');
DmsB  = interp1(azz, Dminus_stream, aB, 'spline');
DTB3D = 0.599*(Deltagro*DpgB + Deltadec*DpdB -fc*(Deltacom + Deltastr*DmsB));

%% There are some cells which have different signs for DTA3D and DTB3D.
%% This can cause the fitting function below to generate complex number!
%% Just for these cells, reset DTA3D. This is acceptable because DTA3D is small
%% in such cases. Best practice is to find a patch with same signs for A & B.
DTA3D = sign(DTA3D).*sign(DTB3D).*DTA3D;

%% Use the fit in appendix of A16.
beta            = 2.8
gamma           = 0.33
alphaYcoeff     = 1/((beta-1)^gamma-(beta-2)^gamma);
powDTA          = (beta-1)^gamma*alphaYcoeff;
powDTB          = (beta-2)^gamma*alphaYcoeff;
alpha3D         = log10(DTB3D./DTA3D)*alphaYcoeff;
coeff_Delta_T3D = 10.^-(   ( (beta-2)^gamma*log10(abs(DTB3D)) -(beta-1)^gamma*log10(abs(DTA3D)) )*alphaYcoeff   );
DT3D_azend      = sign(DTA3D).* 10.^(alpha3D.*(log10(azend)+beta).^gamma) .* coeff_Delta_T3D;
