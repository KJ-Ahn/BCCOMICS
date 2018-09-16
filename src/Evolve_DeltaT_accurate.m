%% We do not use this script, but for accurate DeltaT evolution before
%% decoupling from photons, this is necessary. However, such an accuracy
%% does not seem necessary due to smallness of DeltaT during when the
%% fit in Appendix of A16 is invalid (z>~500).
%% CAMB does not do this calculation, so this script was born.
%% Keeping this for academic reason and possible future improvement.


%% load pre-calcuated radiation density transfer function, with uniform log10(z+1)
load 'CAMB_for_photon/transfer_gamma_table.matbin' '-mat' 'Nout' 'Nk' 'zrecord' 'krecord' 'trk';
global zr_CAMB drk_CAMB log10zr11 dlog10zr1;
krecord   = krecord*h; %% now krecord is in Mpc^-1 unit
zr_CAMB   = zrecord;
azrecord  = 1./(1+zrecord);
%%azrecord_to5 = azrecord(1:lookUP(zrecord,5)); %% gnu octave function lookUP used
Nzrecord     = length(azrecord);
azrecord_to5 = azrecord(1:floor(interp1(zrecord,[1:Nzrecord]',5))); %% mimicing octave lookup function'
zzrecord_to5 = 1./azrecord_to5-1;
Nzrecord_to5 = length(azrecord_to5);
log10zr11    = log10(zrecord(1)+1);
dlog10zr1    = log10(zrecord(1)+1)-log10(zrecord(2)+1);
azbegin      = azrecord_to5(1);
%% azend becomes different from above now.
azend        = azrecord_to5(Nzrecord_to5);

kkmin       = 2*pi/Lbox;
kkmax       = 2*pi/Lcell;
%% ikkmin      = lookUP(ktab,kkmin);
%% ikkmax      = lookUP(ktab,kkmax);
ikkmin      = floor(interp1(ktab,[1:length(ktab)]',kkmin)); %%'
ikkmax      = floor(interp1(ktab,[1:length(ktab)]',kkmax)); %%'
kk4box      = ktab(ikkmin:ikkmax);
Nkk         = length(kk4box);
deltaT_evol = zeros(Nkk,Nzrecord_to5);
%% takes a while to calculate deltaT_evol, so save it after calculation
if (~exist('deltaT_evol.matbin'))
  for ikk = ikkmin:ikkmax
    %% krecord is NOT exactly uniform in log, so need to look up the table.
    %%    indkk = lookUP(krecord, ktab(ikk));
    indkk = floor(interp1(krecord, [1:length(krecord)], ktab(ikk)));
    trk_CAMB = interp1([krecord(indkk) krecord(indkk+1)], [trk(indkk,:); trk(indkk+1,:)], ktab(ikk), 'linear'); %% evolution of deltaTgamma at given k
    Prk_CAMB = As*(ktab(ikk)/k0)^(ns-1) *ktab(ikk) *2*pi^2 *h^3 * trk_CAMB.^2;
    Prk_CAMB = Prk_CAMB * h^-3; %% now in Mpc^3 unit
    drk_CAMB = sqrt(Prk_CAMB).*sign(trk_CAMB);
    Deltagro_kk   = interp1(kktab, Deltagro_k  , ktab(ikk), 'spline');
    Deltadec_kk = interp1(kktab, Deltadec_k, ktab(ikk), 'spline');
    Deltastr_kk   = interp1(kktab, Deltastr_k , ktab(ikk), 'spline');
    x0    = delta_T(ikk);
    options = odeset('RelTol',1e-4,'AbsTol',0.01*min(abs(x0)));
    [azode, dTode] = ode45(@fdDTda,[azbegin,azend],x0,options);
    dT_evol = interp1(azode, dTode, azrecord_to5, 'linear');
    deltaT_evol(ikk-ikkmin+1,:) = dT_evol;
  end
  save('deltaT_evol.matbin', 'azrecord_to5', 'kk4box', 'deltaT_evol', '-v6');
else
  load 'deltaT_evol.matbin' '-mat' 'azrecord_to5' 'kk4box' 'deltaT_evol'
end

%%%% debug--- compare to Naoz & Barkana figure 1
%%%% 
idx1000 = floor(interp1(1./azrecord_to5-1, [1:Nzrecord_to5]', 1000));
idx800  = floor(interp1(1./azrecord_to5-1, [1:Nzrecord_to5]', 800));
idx400  = floor(interp1(1./azrecord_to5-1, [1:Nzrecord_to5]', 400));
idx200  = floor(interp1(1./azrecord_to5-1, [1:Nzrecord_to5]', 200));
idxsparse = [idx1000; idx800; idx400; idx200];
Nsparse   = length(idxsparse)


loglog(kk4box, kk4box.^1.5.*abs(db_evol(:,idx200))/sqrt(2*pi^2), kk4box, kk4box.^1.5.*abs(deltaT_evol(:,idx200))/sqrt(2*pi^2));
print -dpng 'Delta_of_DeltaT_Deltab_z200_approx.png'
close

