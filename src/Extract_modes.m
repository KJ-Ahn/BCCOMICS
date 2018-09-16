%% Extract modes, using algebraic relation (eq. 29 of A16)

%%%%%%%% Prepare to extract 4 modes. ------------------------------------------- begin
%%
%% Extraction of modes using Eq. 29 of Ahn16.
%%
%% Practically z=1000 and z=800 combination works best for k>~ 0.01 Mpc^{-1}.
%% (see Fig. 10 of Ahn16)
%%
%% One may instead use Eq. 7 of Ahn16 algebraically at one single z, but
%% modes extracted this way gives too much error when constructing Delta_{+}
%% and Delta_{-}, so we do NOT do this. 
%%
%% Some thoughts on erros (as in Fig. 10 of Ahn16) -------
%% Fluctuation in radiation components are ignored, but it is not small for
%% k<~k_{eq}~0.01 Mpc^{-1} around recombination. Our objective is to evolve 
%% 4 Mpc volume, so we can safely ignore radiation components. There also can
%% arise potential (e.g. phi and psi in conformal Newtonian gauge) 
%% related terms (see Eq 23a and 43 in MB) in both the continuity equation and 
%% the Poisson equation, but again 4 Mpc is well inside horizon and these terms 
%% are negligible to make Newtonian perturbation theory valid. 
%% Still, the reason why the match is not perfect for large k needs to be
%% understood.
%% -------------------------------------------------------
%%
%% Do this from CAMB outputs under /CAMB_for_mode_finding.
%%
%% z=1000 and 800 are only necessary; other redshifts are to see how well these modes
%% reproduce CAMB-calculated values.
zzz  = load(TFzred);
Nzzz = length(zzz);
azzz = 1./(1+zzz);

%% Check if redshift file correctly describes existing TF files
for izzz=1:Nzzz
  %% see CAMB ReadMe at http://camb.info for meaning of columns
  if ~exist([TFstr1 num2str(zzz(izzz)) TFstr2])
    disp(['File ' TFstr1 num2str(zzz(izzz)) TFstr2 ' does not exist.']);
    disp(['Check if transfer function at z=' num2str(zzz(izzz)) ' exist.']);
    disp(['The redshift description of the transfer file should match']);
    disp(['the format of redshift:' num2str(zzz(izzz))]);
    returnflag = true;
    return;
  end
end

for izzz=1:Nzzz
  %% see CAMB ReadMe at http://camb.info for meaning of columns
  TFF = load([TFstr1 num2str(zzz(izzz)) TFstr2]);

  %% power spectrum
  Pkkc(:,izzz)  = PS_wo_TFtab .* TFF(:,2).^2;  %% Mpc^3 unit, CDM
  Pkkb(:,izzz)  = PS_wo_TFtab .* TFF(:,3).^2;  %% Mpc^3 unit, baryon

  %% perturbation
  Dc(:,izzz)  =  sqrt(Pkkc(:,izzz)) .*sign(TFF(:,2))      ; %% Mpc^(3/2) unit
  Db(:,izzz)  =  sqrt(Pkkb(:,izzz)) .*sign(TFF(:,3))      ; %% Mpc^(3/2) unit

  %% growth factors
  a            = azzz(izzz);
  Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
  Dpg_zz(izzz) = Dpg;
  Dpd_zz(izzz) = Dpd;
  Dms_zz(izzz) = Dms;
end
%%%%%%%% Prepare to extract 4 modes. ------------------------------------------- end

%%%%%%%% Find modes (growing, decaying, compensated, streaming) ---------------- begin
iz1 = lookUP(zzz,1000);
iz2 = lookUP(zzz,800);
if (zzz(iz1)~=1000 || zzz(iz2)~=800)  %% little safeguard
  disp('Designated redshifts for mode extraction not chosen or files nonexixtent.')
  returnflag = true;
  return;
end

Delta_plus_1  = fc*Dc(:,iz1) + fb*Db(:,iz1);
Delta_plus_2  = fc*Dc(:,iz2) + fb*Db(:,iz2);
Delta_minus_1 =    Dc(:,iz1) -    Db(:,iz1);
Delta_minus_2 =    Dc(:,iz2) -    Db(:,iz2);

Dms1 = Dms_zz(iz1);
Dms2 = Dms_zz(iz2);
Dpg1 = Dpg_zz(iz1);
Dpg2 = Dpg_zz(iz2);
Dpd1 = Dpd_zz(iz1);
Dpd2 = Dpd_zz(iz2);

%% Eq. 29 of Ahn16, but generalized for generic z1.
%% One can try different set of zz, but {1000, 800} is found optimal for k>~0.01/Mpc
%% and all redshifts.
Deltagro_k   = (Delta_plus_2*Dpd1 - Delta_plus_1*Dpd2)/(Dpg2*Dpd1 - Dpd2*Dpg1);
Deltadec_k   = (Delta_plus_1 - Deltagro_k*Dpg1)/Dpd1;
Deltastr_k   = (Delta_minus_2 - Delta_minus_1) /(Dms2 - Dms1);
Deltacom_k   = Delta_minus_1 - Deltastr_k*Dms1              ; 
%%%%%%%% Find modes (growing, decaying, compensated, streaming) ---------------- end
