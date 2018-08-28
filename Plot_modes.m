%%%% Script: plots modes
%%

ifig=1;
figure(ifig);
%% plot modes at z=1000
loglog(kktab, Deltagro_k, kktab, Deltadec_k, kktab, -Deltastr_k, kktab, Deltacom_k);
axis([1e-2 10 1e-5 10]);
legend('growing', 'decaying', 'streaming', 'compensated');
title('z=1000');

%% plot modes at z=500
ifig=ifig+1;
figure(ifig);
a = 1/(1+500);
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
loglog(kktab, Deltagro_k*Dpg, kktab, Deltadec_k*Dpd, kktab, -Deltastr_k*Dms, kktab, Deltacom_k);
axis([1e-2 10 1e-5 20]);
legend('growing', 'decaying', 'streaming', 'compensated');
title('z=500');

%% plot modes at z=200
ifig=ifig+1;
figure(ifig);
a = 1/(1+200);
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
loglog(kktab, Deltagro_k*Dpg, kktab, Deltadec_k*Dpd, kktab, -Deltastr_k*Dms, kktab, Deltacom_k);
axis([1e-2 10 1e-5 30]);
legend('growing', 'decaying', 'streaming', 'compensated');
title('z=200');

%% Compare CAMB output to fluctuations constructed from extracted modes
ifig=ifig+1;
figure(ifig);
izzz        = lookUP(zzz, 500);
Dplus_CAMB  = fc*Dc(:,izzz) + fb*Db(:,izzz);
Dminus_CAMB =    Dc(:,izzz) -    Db(:,izzz);
a           = 1/(1+zzz(izzz));
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
Dplus_cons  = Deltagro_k*Dpg + Deltadec_k*Dpd;
Dminus_cons = Deltacom_k     + Deltastr_k*Dms;
semilogx(kktab, Dplus_cons./Dplus_CAMB-1, kktab, Dminus_cons./Dminus_CAMB-1);
axis([1e-3 100 -0.1 0.05])
legend('\Delta_{+}','\Delta_{-}')
title('z=500; fractional difference');

%% highest k difference expected because of (ignored in Eq. 5 of Ahn16) pressure term.
ifig=ifig+1;
figure(ifig);
izzz        = lookUP(zzz, 50);
Dplus_CAMB  = fc*Dc(:,izzz) + fb*Db(:,izzz);
Dminus_CAMB =    Dc(:,izzz) -    Db(:,izzz);
a           = 1/(1+zzz(izzz));
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
Dplus_cons  = Deltagro_k*Dpg + Deltadec_k*Dpd;
Dminus_cons = Deltacom_k     + Deltastr_k*Dms;
semilogx(kktab, Dplus_cons./Dplus_CAMB-1, kktab, Dminus_cons./Dminus_CAMB-1);
axis([1e-3 100 -0.1 0.05])
legend('\Delta_{+}','\Delta_{-}')
title('z=50; fractional difference');

%% highest k difference expected because of (ignored in Eq. 5 of Ahn16) pressure term.
ifig=ifig+1;
figure(ifig);
izzz        = lookUP(zzz, 5);
Dplus_CAMB  = fc*Dc(:,izzz) + fb*Db(:,izzz);
Dminus_CAMB =    Dc(:,izzz) -    Db(:,izzz);
a           = 1/(1+zzz(izzz));
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
Dplus_cons  = Deltagro_k*Dpg + Deltadec_k*Dpd;
Dminus_cons = Deltacom_k     + Deltastr_k*Dms;
semilogx(kktab, Dplus_cons./Dplus_CAMB-1, kktab, Dminus_cons./Dminus_CAMB-1);
axis([1e-3 100 -0.1 0.05])
legend('\Delta_{+}','\Delta_{-}')
title('z=5; fractional difference');

%% highest k difference expected because of (ignored in Eq. 5 of Ahn16) pressure term.
ifig=ifig+1;
figure(ifig);
izzz        = lookUP(zzz, 0);
Dplus_CAMB  = fc*Dc(:,izzz) + fb*Db(:,izzz);
Dminus_CAMB =    Dc(:,izzz) -    Db(:,izzz);
a           = 1/(1+zzz(izzz));
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
Dplus_cons  = Deltagro_k*Dpg + Deltadec_k*Dpd;
Dminus_cons = Deltacom_k     + Deltastr_k*Dms;
semilogx(kktab, Dplus_cons./Dplus_CAMB-1, kktab, Dminus_cons./Dminus_CAMB-1);
axis([1e-3 100 -0.1 0.05])
legend('\Delta_{+}','\Delta_{-}')
title('z=0; fractional difference');
