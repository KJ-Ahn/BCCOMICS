%% BCCOMICS_CR_SETUP: Sets up patch values at recombination using Constrained Realization (CR).

clear;
more off; %% enables to see progress
returnflag = false; %% main program need to stop when script stops.

%% Start recording log
diary on;

%% Detect which is running: octave or matlab?
if (exist('OCTAVE_VERSION','builtin'))
  matlabflag=false;
  disp('----------------run on OCTAVE----------------');
else
  matlabflag=true;
  disp('----------------run on MATLAB----------------');
end

disp('----------------Initializing----------------');

%% Read in essential parameters
run('params.m');  %%==== script ==================
%% Define some box-related quantities
box_init;  %%==== script ==================
if (mod(Ncell,2)==0)
  disp('Choose an odd number to make a patch 4 Mpc in size');
  clear; %% clearing workspace (memory)
  return;
end
if (Lcell ~= 4)
  disp('Make Lcell as close as to 4 Mpc; otherwise DeltaT will gain some error.');
  disp('Will proceed anyway, but you have been warned...');
  disp(' --- Hit any key to proceed, or Ctrl+C to stop --- ');
  pause;
end

%% Some old versions of gnu octave has buggy ifftshift routine, so for
%% octave version older than 4.0.1, just use working one under the provided
%% directory. In case ode45 is not available, use provided ODE package. 
%% All this can be avoided by upgrading to most recent octave version.
if ~matlabflag
  if compare_versions(OCTAVE_VERSION,'4.0.1','<')
    %% Messages "warning: function * shadows ..." should be welcomed.
    addpath([pkgdir '/mfiles_for_octave']); 
  end
  if ~exist('ode45')
    addpath([pkgdir '/odepkg-0.8.5']);
  end
end

%% Create directory to dump outputs
if ~exist(setupdir)
  mkdir(setupdir);
end

%% constants
global mH kb MpcMyr_2_kms;
%% Read in constants in cgs unit and conversion factors.
Consts_Conversions;  %%==== script ==================

%% cosmological parameters
global H0 Om0 Omr0 TCMB0 OmLambda0;
global tgamma;
%% read in cosmological parameters for background LambdaCDM universe
%% -- CAREFUL: Numerical values need to match CAMB input !!!!!!!!!!
run(Cosmology);  %%==== script ==================

%% Working at only very high z, so below is OK for now.
global fb fc;
fb = ombh2/(ombh2+omch2); %% baryon/matter fraction
fc = omch2/(ombh2+omch2); %% CDM/matter fraction

%% Using fit by TH (Eq. 2) for global baryon temperature.
%% Do NOT change zi below.
%% Also the initial transfer function is loaded here.
global ai aa1 aa2;
zi      = 1000;  %% our choice for beginning redshift (soon after recombination)
ai      = 1/(1+zi);
Hzi     = H0*sqrt(Om0*(1+zi)^3 + Omr0*(1+zi)^4); %% initil Hubble in Myr^-1 unit
aa1     = 1/119; %% aa1 & aa2 under Eq. 2 in TH
aa2     = 1/115;
Tbzi    = TCMB0/ai /(1+ai/aa1/(1+(aa2/ai)^1.5)); %% baryon temperature fit in NB
cszi    = sqrt((5/3)*kb*Tbzi/(1.22*mH)) * 1e-5;  %% sound speed in km/s
Tgammai = TCMB0*(1+zi);  %% CMB temperature at z=1000

%% Fluctuations and power spectra at zi ----------------------- begin
TF_zi = load([TFstr1 num2str(zi) TFstr2]); %% transfer function at zi
kktab = TF_zi(:,1)*h;  %% k, in Mpc^-1 unit

%% Primordial power spectrum: See IV.A in CAMB.pdf from http://cosmologist.info/notes
lnPstab     = log(As)+(ns-1)*log(kktab/k0)+nrun/2*(log(kktab/k0)).^2+nrunrun/6*(log(kktab/k0)).^3;
%% powe spectrum without TF^2, where TF is the CAMB transfer function output
%% Refer to Transfer_GetMatterPowerData subroutine in CAMB
PS_wo_TFtab_ = exp(lnPstab) .* kktab *2*pi^2 * h^3; %% if TF^2 multiplied, in h^-3 Mpc^3 unit
PS_wo_TFtab  = PS_wo_TFtab_ * h^-3;  %% if TF^2 multiplied, in Mpc^3 unit

Pkc_zi   = PS_wo_TFtab .* TF_zi(:,2).^2;  %% Mpc^3 unit, CDM
Pkb_zi   = PS_wo_TFtab .* TF_zi(:,3).^2;  %% Mpc^3 unit, baryon
Pkr_zi   = PS_wo_TFtab .* TF_zi(:,4).^2;  %% Mpc^3 unit, radiation
PkTHc_zi = PS_wo_TFtab .* TF_zi(:,11).^2;  %% Mpc^3 unit, CDM vel divergence
PkTHb_zi = PS_wo_TFtab .* TF_zi(:,12).^2;  %% Mpc^3 unit, baryon vel divergence
PkVcb_zi = PS_wo_TFtab .* TF_zi(:,13).^2;  %% Mpc^3 unit, Vc-Vb

%% perturbation -- see CAMB Readme for meaning of columns
Dc_zi  =  sqrt(Pkc_zi)   .*sign(TF_zi(:,2)); %% Mpc^(3/2) unit
Db_zi  =  sqrt(Pkb_zi)   .*sign(TF_zi(:,3)); %% Mpc^(3/2) unit
Dr_zi  =  sqrt(Pkr_zi)   .*sign(TF_zi(:,4)); %% Mpc^(3/2) unit
THc_zi = -sqrt(PkTHc_zi) .*sign(TF_zi(:,11))*Hzi; %% Mpc^(3/2) Myr^-1 unit
THb_zi = -sqrt(PkTHb_zi) .*sign(TF_zi(:,12))*Hzi; %% Mpc^(3/2) Myr^-1 unit
Vcb_zi = -sqrt(PkVcb_zi) .*sign(TF_zi(:,13))*c_inkms/MpcMyr_2_kms; %% Mpc^(3/2) Mpc Myr^-1 unit

%% sanity check of sign: Try plots for confirmation if wanted...
Vc_zi  = -sqrt(PkTHc_zi) .*sign(TF_zi(:,11))*ai*Hzi./kktab;
Vb_zi  = -sqrt(PkTHb_zi) .*sign(TF_zi(:,12))*ai*Hzi./kktab;
% semilogx(kktab, THc_zi); %% should be negative, because THc = -dDc/dt
% semilogx(kktab, (Vc_zi-Vb_zi)./Vcb_zi); %% should be 1, NOT -1.
%% Fluctuations and power spectra at zi ----------------------- end

%% Use recfast output for z(redshift)-xe(global ionized fraction) table.
%% Table can be non-recfast as long as you trust it.
%% CAREFUL: for efficient interpolation at any given redshift, the redshift
%% interval MUST be uniform.
global zrecf xerecf dzrecf zrecf1; %% do not bother naming convention
zxe    = load(zxestr); %% Well recfast is preferred and thus the name of variables...
zrecf  = zxe(:,1);
xerecf = zxe(:,2);
xei    = interp1(zrecf, xerecf, zi); %% xe(zi)
dzrecf = zrecf(1)-zrecf(2);
zrecf1 = zrecf(1);

%% a bit of safeguard for non-uniform-z recfast table
if (dzrecf ~= zrecf(9)-zrecf(10))
  disp('Recfast output is not uniform in z. Quitting.');
  clear;
  return;
end
%% table should be in descending order in z
if (zrecf(1) < zrecf(2))
  disp('z-xe table should have decreasing z');
  clear;
  return;
end

%% Temporal evolution of growing, decaying, and streaming modes ---------------- begin
%% Radiation components (photon + neutrino) make these modes NOT follow the simple
%% power laws (\propto a, a^-1.5, a^-0.5 respectively for density). Therefore,
%% numerical integration should be done to find the correct mode evolution &
%% mode extraction.
global Dplus_grow Dplus_decay Dminus_stream dDplus_grow_da dDplus_decay_da dDminus_stream_da azz log10az_min dlog10az;
Get_growth;  %%==== script ==================
%% Temporal evolution of growing, decaying, and streaming modes ---------------- end

%%%%%%%% mode extraction and, if wanted, some plotting ----------------------- begin
Extract_modes;  %%==== script ==================
if returnflag %% inherit scriptwise return
  clear;
  return;
end
  %% --- plot ---
if plotflag
  Plot_modes;  %%==== script ==================
end
%% dump modes
fout  = fopen([setupdir '/k_gro_dec_com_str.dat'],'w');
mdata = [kktab Deltagro_k Deltadec_k Deltacom_k Deltastr_k];
fprintf(fout,'%14.7e %14.7e %14.7e %14.7e %14.7e\n', mdata'); %%'
fclose(fout);
%%%%%%%% mode extraction and, if wanted, some plotting ------------------------- end

%% At z=1000, get baryon temperature fluctuation DT_zi
Initialize_DeltaT;  %%==== script ==================

%%%%%% Get 3D spatial fluctuatons in the big box ------------------------------ begin
%% 2 Gaussian random number sets into a complex number field.
%% Generate AND use file only when there does NOT exist the random seed file.
fgaussstr = [setupdir '/gaussseed.matbin'];
if (~exist(fgaussstr))
  disp(['------gaussian seed ' fgaussstr ' being generated------']);
  Nreallization = Ncell^3;  %% need to be >= Ncell^3
  gauss1 = normrnd(0,1,[Nreallization,1]);
  gauss2 = normrnd(0,1,[Nreallization,1]);
  gauss  = gauss1+i*gauss2;
  %% Save gaussian ramdom seed, in complex format.
  %% For compatibility with older octave versions, matlab binary should be in v6.
  if (matlabflag)
    save(fgaussstr, 'gauss', '-v6');
  else
    save('-mat-binary', fgaussstr, 'gauss'); %% -mat-binary = -v6 in octave
  end
else %% load preexisting seed
  disp(['------preexising gaussian seed ' fgaussstr ' being read in------']);
  if matlabflag
    load(fgaussstr, '-mat', 'gauss')  %% matlab knows about binary version
  else
    load('-mat-binary', fgaussstr, 'gauss')
  end
end
%% 3D gaussian random number (G1+iG2)
gauss3D = reshape(gauss(1:Nmode^3),Nmode,Nmode,Nmode);

disp('----------------Generating real-space 3D fields----------------');
%% Random-number-seeding and FFTing to generate 3D fields of patches at z=zi=1000
Get_patches_3D_zi;  %%==== script ==================
%%%%%% Get 3D spatial fluctuatons in the big box ------------------------------ end

%%%%%% Constrained Realization Pipeline --------------------------------------- begin
%% Calculate background statistics from the unconstrained realization
disp('----------------Calculating Base Statistics----------------');
Delta_m = fc * Delta_c + fb * Delta_b;
stdDm   = std(Delta_m(:));
rmsVcb  = sqrt(mean(V_cb_1(:).^2 + V_cb_2(:).^2 + V_cb_3(:).^2));

%% Phase 1: Set Target Constraints
Set_constraints;  %%==== script ==================

%% Phase 2: Calculate Covariance Matrix for CR
Calculate_Covariance;  %%==== script ==================

%% Phase 3: Apply Constrained Realization Filter to Fourier Space
Apply_CR_filter;  %%==== script ==================
%%%%%% Constrained Realization Pipeline --------------------------------------- end

azbegin = ai;    %% z=1000
azend   = 1/(1+zzend);

%% By using modes generate 3D fields of patches at z=zzend.
%% Only for DeltaT, fitting formula is used (Get_DeltaT_fit used inside the script)
global beta gamma;
global signDT alpha coeff_Delta_T;
Get_patches_3D_zend;  %%==== script ==================

disp('----- Recording Background Statistics (stats_zi.dat, stats_zend.dat) -----');
stdDc   = std(Delta_c(:));
stdTc   = std(Theta_c(:));
stdDb   = std(Delta_b(:));
stdTb   = std(Theta_b(:));
rmsVc   = sqrt(mean(V_c_1(:).^2+V_c_2(:).^2+V_c_3(:).^2));
rmsVb   = sqrt(mean(V_b_1(:).^2+V_b_2(:).^2+V_b_3(:).^2));
Vcbp    = sqrt(2/3)*rmsVcb;
stdT    = std(Delta_T(:));
datstat_zi = [stdDc stdTc rmsVc*MpcMyr_2_kms stdDb stdTb rmsVb*MpcMyr_2_kms rmsVcb*MpcMyr_2_kms Vcbp*MpcMyr_2_kms stdT stdDm];
fout    = fopen([setupdir '/stats_zi.dat'],'w');
fprintf(fout,'@ zi:  stdDc stdThc rmsVc stdDb stdThb rmsVb rmsVcb Vcbp stdT stdDm\n');
fprintf(fout,'Units: None  Myr^-1 km/s  None  Myr^-1 km/s  km/s   km/s None None\n');
fprintf(fout,'%e %e %e %e %e %e %e %e %e %e\n', datstat_zi);
fclose(fout);

sDc_azend    = std(Dc3D_azend(:));
sTc_azend    = std(THc3D_azend(:));
Dm3D_azend   = fc * Dc3D_azend + fb * Db3D_azend;
sDm_azend    = std(Dm3D_azend(:));
sDb_azend    = std(Db3D_azend(:));
sTb_azend    = std(THb3D_azend(:));
Vcb_azend    = sqrt(V_cb_1_azend(:).^2+V_cb_2_azend(:).^2+V_cb_3_azend(:).^2);
rmsVcb_azend = sqrt(mean(Vcb_azend.^2));
Vcbp_azend   = sqrt(2/3)*rmsVcb_azend;
sT_azend     = std(DT3D_azend(:));
datstat_zend = [sDc_azend sTc_azend sDb_azend sTb_azend rmsVcb_azend*MpcMyr_2_kms Vcbp_azend*MpcMyr_2_kms sT_azend sDm_azend];
fout    = fopen([setupdir '/stats_zend.dat'],'w');
fprintf(fout,'@ zend:  stdDc  stdThc  stdDb  stdThb  rmsVcb  Vcbp  stdT sDm_azend\n');
fprintf(fout,'Units:   None   Myr^-1  None   Myr^-1  km/s    km/s  None None\n');
fprintf(fout,'%e %e %e %e %e %e %e %e\n', datstat_zend);
fclose(fout);
%% =======================================================================

%%%%%% Data Dumping for bccomics_CR.m -------------------------------------------- begin
disp('----- Saving CR-filtered 3D fields at z=1000 -----');
if matlabflag
  save([setupdir '/Dc3D.matbin'],   'Delta_c', '-v6'); 
  save([setupdir '/Db3D.matbin'],   'Delta_b', '-v6'); 
  save([setupdir '/THc3D.matbin'],  'Theta_c', '-v6'); 
  save([setupdir '/THb3D.matbin'],  'Theta_b', '-v6'); 
  save([setupdir '/DT.matbin'],     'Delta_T', '-v6');
  save([setupdir '/V_cb_1.matbin'], 'V_cb_1',  '-v6'); 
  save([setupdir '/V_cb_2.matbin'], 'V_cb_2',  '-v6'); 
  save([setupdir '/V_cb_3.matbin'], 'V_cb_3',  '-v6'); 
else
  save('-mat-binary', [setupdir '/Dc3D.matbin'],   'Delta_c'); 
  save('-mat-binary', [setupdir '/Db3D.matbin'],   'Delta_b'); 
  save('-mat-binary', [setupdir '/THc3D.matbin'],  'Theta_c'); 
  save('-mat-binary', [setupdir '/THb3D.matbin'],  'Theta_b'); 
  save('-mat-binary', [setupdir '/DT.matbin'],     'Delta_T');
  save('-mat-binary', [setupdir '/V_cb_1.matbin'], 'V_cb_1'); 
  save('-mat-binary', [setupdir '/V_cb_2.matbin'], 'V_cb_2'); 
  save('-mat-binary', [setupdir '/V_cb_3.matbin'], 'V_cb_3'); 
end

disp('----- Recording target patch data -----');
%% Record the single target coordinate (Append mode)
fout = fopen([setupdir '/icc.dat'], 'a');
fprintf(fout, '%i %i %i\n', icc');
fclose(fout);

%% Data to dump MUST remain as Delta_c / Delta_b for subsequent steps!
daticc_zi(1) = Delta_c(icc(1),icc(2),icc(3));
daticc_zi(2) = Delta_b(icc(1),icc(2),icc(3));
daticc_zi(3) = Theta_c(icc(1),icc(2),icc(3));
daticc_zi(4) = Theta_b(icc(1),icc(2),icc(3));
daticc_zi(5) = V_cb_1 (icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zi(6) = V_cb_2 (icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zi(7) = V_cb_3 (icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zi(8) = norm([V_cb_1(icc(1),icc(2),icc(3)) V_cb_2(icc(1),icc(2),icc(3)) V_cb_3(icc(1),icc(2),icc(3))])*MpcMyr_2_kms;
daticc_zi(9) = Delta_T(icc(1),icc(2),icc(3));
fout = fopen([setupdir '/zi_icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat'],'a');
fprintf(fout,'%i %i %i %e %e %e %e %e %e %e %e %e\n',[icc daticc_zi]');
fclose(fout);

%% Record zend data (utilizing variables updated by Get_patches_3D_zend)
daticc_zend(1) = Dc3D_azend  (icc(1),icc(2),icc(3));
daticc_zend(2) = Db3D_azend  (icc(1),icc(2),icc(3));
daticc_zend(3) = THc3D_azend (icc(1),icc(2),icc(3));
daticc_zend(4) = THb3D_azend (icc(1),icc(2),icc(3));
daticc_zend(5) = V_cb_1_azend(icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zend(6) = V_cb_2_azend(icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zend(7) = V_cb_3_azend(icc(1),icc(2),icc(3))*MpcMyr_2_kms;
daticc_zend(8) = norm([V_cb_1_azend(icc(1),icc(2),icc(3)) V_cb_2_azend(icc(1),icc(2),icc(3)) V_cb_3_azend(icc(1),icc(2),icc(3))])*MpcMyr_2_kms;
daticc_zend(9) = DT3D_azend  (icc(1),icc(2),icc(3));
fout = fopen([setupdir '/zend_icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat'],'a');
fprintf(fout,'%i %i %i %e %e %e %e %e %e %e %e %e\n',[icc daticc_zend]');
fclose(fout);
%%%%%% Data Dumping for bccomics_CR.m -------------------------------------------- end

%% master equation for high k modes:
%% *_p are the 4 modes at a chosen patch
global ksample costh Deltagro_p Deltadec_p Deltacom_p Deltastr_p;
global Thc_i Thb_i rV_i;

%% Assigning initial target variables at chosen patch ---------- begin
ic = icc(1);
jc = icc(2);
kc = icc(3);

Deltagro_p = Deltagro(ic,jc,kc);
Deltadec_p = Deltadec(ic,jc,kc);
Deltacom_p = Deltacom(ic,jc,kc);
Deltastr_p = Deltastr(ic,jc,kc);

Thc_i = Theta_c(ic,jc,kc);
Thb_i = Theta_b(ic,jc,kc);
rV_i  = norm([V_cb_1(ic,jc,kc), V_cb_2(ic,jc,kc), V_cb_3(ic,jc,kc)]);
%% Assigning initial target variables at chosen patch ---------- end

%% for mu(=cosine of angle between Vcb and k) loop
dmu = 0.05;
mu  = 0:dmu:1; %% Use symmetry of P(k,mu) about mu=0 to save calculation time.
Nmu = length(mu);

if matlabflag
  save([setupdir '/mu.dat'],'mu','-ascii');
else
  save('-ascii',[setupdir '/mu.dat'],'mu');
end

disp('----------------Integrating----------------');

Integrate_evolODE;  %%==== script ==================

disp('*********** bccomics_CR_setup successfully ended ************');
diary off;
movefile('diary', 'bccomics_CR_setup.log');
