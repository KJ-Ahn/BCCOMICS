%% Integrate evolution equation (eq. 11 of A16) and dump its output

%% Loop over different k values ----------------------------------------begin
ksampletab = [[1:28]'; 30; 33; 40; 45; 50; 56; 60; 75; 80; 90; 100; 110; 125; 135; 155; 170; 185; 200; 225; 260; 290; 335; 360; 400; 430; 460; 500; 530; 560; 600; 640; 680; 720; 760; 800; 840; 880; 930; 980; 1000; 1100; 1200; 1300; 1400; 1550; 1600; 1700; 1800; 1900; 2000; 3500; 5600; 9300; 13500];
Nsample    = length(ksampletab)

%% directory for main output (fluctuation)
TFdir = [outputdir '/deltas'];
if ~exist(TFdir)
  mkdir(TFdir);
end
if matlabflag
  save([TFdir '/ksample.dat'], 'ksampletab', '-ascii');
else
  save('-ascii', [TFdir '/ksample.dat'], 'ksampletab');
end

TFTHdir = [outputdir '/deltas_TH'];
if THflag
  if ~exist(TFTHdir)
    mkdir(TFTHdir);
  end
  if matlabflag
    save([TFTHdir '/ksample.dat'], 'ksampletab', '-ascii');
  else
    save('-ascii', [TFTHdir '/ksample.dat'], 'ksampletab');
  end
end

%% high-k mode fluctuations at z=zi.
delta_c_k = interp1(kktab, Dc_zi,  ksampletab, 'spline');
theta_c_k = interp1(kktab, THc_zi, ksampletab, 'spline');
delta_b_k = interp1(kktab, Db_zi,  ksampletab, 'spline');
theta_b_k = interp1(kktab, THb_zi, ksampletab, 'spline');
delta_T_k = interp1(kktab, DT_zi,  ksampletab, 'spline');

%% Time integration interval for ode45
azbegin = ai;    %% z=1000
azend   = 1/(1+zzend);
zzbegin = 1/azbegin-1;
zzend   = 1/azend  -1;

if matlabflag
  save([outputdir '/zz.dat'], 'zzbegin','zzend', '-ascii');
else
  save('-ascii', [outputdir '/zz.dat'], 'zzbegin','zzend');
end

ic = icc1;
jc = icc2;
kc = icc3;
for isample=1:Nsample  %% this loop takes a long while!!
  disp([num2str(isample) 'th wavenumber out of ' num2str(Nsample) ' is being handled.']);
  ksample    = ksampletab(isample);

  %% Integrate evolution equation by TH, just for comparison.
  x0 = [delta_c_k(isample); theta_c_k(isample); delta_b_k(isample); theta_b_k(isample); delta_c_k(isample); theta_c_k(isample); delta_b_k(isample); theta_b_k(isample); delta_T_k(isample); delta_T_k(isample)]/sqrt(2);

  %% Main output: fluctuation at zzend as a function of k.
  stroutD    = [TFdir '/Deltas_1Dmu_k' num2str(ksample)];
  stroutD    = [stroutD   '_ic' num2str(ic) '_jc' num2str(jc) '_kc' num2str(kc) '-muhalf.matbin'];

  strTHoutD  = [TFTHdir '/Deltas_TH_1Dmu_k' num2str(ksample)];
  strTHoutD  = [strTHoutD '_ic' num2str(ic) '_jc' num2str(jc) '_kc' num2str(kc) '-muhalf.matbin'];

  if exist(stroutD)
    disp([stroutD ' exists. Possibly for all other ks''']);
    return;
  end
  if (~exist(stroutD) || OWRTflag)
    Deltagro_kk   = Deltagro(ic, jc, kc);
    Deltadec_kk   = Deltadec(ic, jc, kc);
    Deltacom_kk   = Deltacom(ic, jc, kc);
    Deltastr_kk   = Deltastr(ic, jc, kc);
    Thc_i         = Theta_c (ic, jc, kc);
    Thb_i         = Theta_b (ic, jc, kc);
    signDT        = sign(DTA3D(ic, jc, kc));
    alpha         = log10(DTB3D(ic,jc,kc)/DTA3D(ic,jc,kc))*alphaYcoeff;
    coeff_Delta_T = abs(DTA3D(ic,jc,kc))^powDTA /abs(DTB3D(ic,jc,kc))^powDTB; %% 10^-Y, Y defined in eq. 30 of A16
    rV_i          = Vcb     (ic, jc, kc);

    %% These are complex values.
    %% for a given patch(icc).
    deltasc   = zeros(Nmu,1); 
    deltasb   = zeros(Nmu,1); 
    deltasThc = zeros(Nmu,1); 
    deltasThb = zeros(Nmu,1);
    deltasT   = zeros(Nmu,1);

    if THflag
      deltasc_TH    = zeros(Nmu,1);
      deltasb_TH    = zeros(Nmu,1); 
      deltasThc_TH  = zeros(Nmu,1); 
      deltasThb_TH  = zeros(Nmu,1);
      deltasT_TH    = zeros(Nmu,1);
    end

    for imu=1:Nmu
      costh = mu(imu);

      %% A16 ================
      options = odeset('RelTol',1e-4,'AbsTol',0.01*min(abs(x0)));
      [azode, deltaode] = ode45(@f,[azbegin,azend],x0,options);
      deltasAhn  = interp1(azode, deltaode, az1, 'linear', 'extrap');
        
      deltasc  (imu) = deltasAhn(Nzz1,1)+i*deltasAhn(Nzz1,5);
      deltasb  (imu) = deltasAhn(Nzz1,3)+i*deltasAhn(Nzz1,7);
      deltasThc(imu) = deltasAhn(Nzz1,2)+i*deltasAhn(Nzz1,6);
      deltasThb(imu) = deltasAhn(Nzz1,4)+i*deltasAhn(Nzz1,8);
      deltasT  (imu) = deltasAhn(Nzz1,9)+i*deltasAhn(Nzz1,10);

      if (THflag && (~exist(strTHoutD) || OWRTflag))
	%% TH ===================
	options = odeset('RelTol',1e-4,'AbsTol',0.01*min(abs(x0)));
	[azode, deltaode] = ode45(@fTH,[azbegin,azend],x0,options);
	deltasTH  = interp1(azode, deltaode, az1, 'linear', 'extrap');
        
	deltasc_TH  (imu) = deltasTH(Nzz1,1)+i*deltasTH(Nzz1,5);
	deltasb_TH  (imu) = deltasTH(Nzz1,3)+i*deltasTH(Nzz1,7);
	deltasThc_TH(imu) = deltasTH(Nzz1,2)+i*deltasTH(Nzz1,6);
	deltasThb_TH(imu) = deltasTH(Nzz1,4)+i*deltasTH(Nzz1,8);
	deltasT_TH  (imu) = deltasTH(Nzz1,9)+i*deltasTH(Nzz1,10);
      end
    end
    if matlabflag
      save(strTHoutD, 'deltasc_TH', 'deltasb_TH', 'deltasThc_TH', 'deltasThb_TH', 'deltasT_TH', '-v6');
      save(stroutD,   'deltasc',    'deltasb',    'deltasThc',    'deltasThb',    'deltasT',    '-v6');
    else
      save('-mat-binary', strTHoutD, 'deltasc_TH', 'deltasb_TH', 'deltasThc_TH', 'deltasThb_TH', 'deltasT_TH');
      save('-mat-binary', stroutD,   'deltasc',    'deltasb',    'deltasThc',    'deltasThb',    'deltasT'   );
    end
  end  %% if-end
end  %% isample loop end

