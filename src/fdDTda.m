%% Part of BCCOMICS.
%%
%% Author: Kyungjin Ahn
%%
%% References - 
%%   TH10: Tseliakhovich & Hirata 2010, PRD 82, 083520
%%   A16: Ahn 2016, ApJ 830:68
%%   AS18: Ahn & Smith 2018, arXiv:1807.04063
%%
%% Rate equation for a Delta_T of a given patch. x = Delta_T. (Eq. 4 of A16)
%% In /da instead of /dt, where a is the scale factor.
%%
%% This equation does NOT drop delta_T_gamma as in Eq. 4 of A16. So it is
%% quite accurate. We extract delta_T_gamma(a) for any scale factor a from
%% CAMB. With drk = 4*delta_T_gamma(a), drk_CAMB is the fluctuation in
%% photon density, calculated by CAMB.
%%


function dxda = fdDTda(a,x)
  global H0 Om0 Omr0 OmLambda0 TCMB0 ai aa1 aa2 fb fc;
  global zrecf xerecf dzrecf zrecf1 tgamma;
  global zr_CAMB drk_CAMB log10zr11 dlog10zr1;
  global dDplus_grow_da dDplus_decay_da dDminus_stream_da azz log10az_min dlog10az;
  global Deltagro_kk Deltadec_kk Deltastr_kk;
  global log10az_min dlog10az

  Hz      = H0*sqrt(Om0*a^-3+Omr0*a^-4 + OmLambda0);
  aH      = a*Hz;

  Tbz     = TCMB0/a /(1+a/aa1/(1+(aa2/a)^1.5));
  Tgamma  = TCMB0/a;

  %% interpolation of deltaTgamma(z) from a table, done efficiently without calling
  %% interp1 with the whole array.
  zz       = 1/a -1;
  Dlog10zr = log10zr11-log10(zz+1);
  indzr1   = floor(Dlog10zr/dlog10zr1)+1;
  indzr2   = indzr1+1;
  zr1      = zr_CAMB(indzr1);
  zr2      = zr_CAMB(indzr2);
  drk1     = drk_CAMB(indzr1);
  drk2     = drk_CAMB(indzr2);
  drk      = interp1([zr1;zr2],[drk1;drk2],zz,'linear','extrap');

  %% interpolation of xe(z) from a table, done efficiently without calling
  %% interp1 with the whole array of zrecf & xerecf.
  ind1    = floor((zrecf1-zz)/dzrecf)+1;
  ind2    = ind1+1;
  zz1     = zrecf(ind1);
  zz2     = zrecf(ind2);
  xe1     = xerecf(ind1);
  xe2     = xerecf(ind2);
  xez     = interp1([zz1;zz2],[xe1;xe2],zz,'linear','extrap');

  %% Patch baryon temperature fluctuation: A16, Eq. 4, 
  %%                                       with drk = delta_gamma = 4*delta_T_gamma
  %% 2/3/a*(...*(a/ai)^-0.5) is 2/3\frac{\partial\Delta_b}{\partial t}/(aH)

  ind1    = floor((log10(a)-log10az_min)/dlog10az)+1;
  ind2    = ind1+1;
  az1     = azz(ind1);
  az2     = azz(ind2);
  dDpg_da = interp1([az1;az2],[dDplus_grow_da(ind1)   ;dDplus_grow_da(ind2)   ],a,'linear','extrap');
  dDpd_da = interp1([az1;az2],[dDplus_decay_da(ind1)  ;dDplus_decay_da(ind2)  ],a,'linear','extrap');
  dDms_da = interp1([az1;az2],[dDminus_stream_da(ind1);dDminus_stream_da(ind2)],a,'linear','extrap');

  dDeltab_da = Deltagro_kk*dDpg_da + Deltadec_kk*dDpd_da - fc*Deltastr_kk;
  dxda       = 2/3*dDeltab_da +xez/tgamma/a^4*(drk*(5/4*Tgamma/Tbz-1)-Tgamma/Tbz*x)/aH;

end

