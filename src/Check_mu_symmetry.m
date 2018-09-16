%% --- This script is NOT used ------------------------------------------
%% Checks if P(k,mu) is symmetric in mu, such that P(k,mu)=P(k,-mu).
%% This fact enables using only mu>=0.

mufull = -1:dmu:1;
Nmufull = length(mufull);

%% just to check symmetry-------------------------------------------------------bebin
deltasc_Ahn_full   = zeros(1, Nmufull); 
deltasb_Ahn_full   = zeros(1, Nmufull); 
deltasThc_Ahn_full = zeros(1, Nmufull); 
deltasThb_Ahn_full = zeros(1, Nmufull);
deltasT_Ahn_full   = zeros(1, Nmufull);

ksample = ksampletab(1);

%% Apply some arbitrary phase, but with coherence over variables.
%% Recommended phase=pi/4
%% Whatever the phase, the power (abs(delta)^2) should look identical, and it is so indeed.
phase=pi/4;
cosph=cos(phase);sinph=sin(phase);
isample=4
x0 = [cosph*[delta_c_k(isample); theta_c_k(isample); delta_b_k(isample); theta_b_k(isample)]; sinph*[delta_c_k(isample); theta_c_k(isample); delta_b_k(isample); theta_b_k(isample)]; cosph*delta_T_k(isample); sinph*delta_T_k(isample)];

ic=1;jc=1;kc=1;
Deltagro_kk   = Dgro3D (ic, jc, kc);
Deltadec_kk   = Ddec3D (ic, jc, kc);
Deltacom_kk   = Dcom3D (ic, jc, kc);
Deltastr_kk   = Dstr3D (ic, jc, kc);
Thc_i         = Theta_c(ic, jc, kc);
Thb_i         = Theta_b(ic, jc, kc);
signDT        = sign(DTA3D(ic, jc, kc));
alpha         = log10(DTB3D(ic,jc,kc)/DTA3D(ic,jc,kc))*alphacoeff;
coeff_Delta_T = abs(DTA3D(ic,jc,kc))^4.259097 /abs(DTB3D(ic,jc,kc))^3.259097;
rV_i          = relV_MpcMyr(sub2ind(size(V_cb_1),ic,jc,kc));

for imu=1:Nmufull
  costh = mufull(imu);
  
  %% Ahn =============================================================================
  options = odeset('RelTol',1e-4,'AbsTol',0.0001*min(abs(x0)));
  [azode, deltaode] = ode45(@f,[azbegin,azend],x0,options);
  deltasAhnfull  = interp1(azode, deltaode, az1, 'linear');
  
  deltasc_Ahn_full(imu)   = deltasAhnfull(Nzz1,1)+i*deltasAhnfull(Nzz1,5);
  deltasb_Ahn_full(imu)   = deltasAhnfull(Nzz1,3)+i*deltasAhnfull(Nzz1,7);
  deltasThc_Ahn_full(imu) = deltasAhnfull(Nzz1,2)+i*deltasAhnfull(Nzz1,6);
  deltasThb_Ahn_full(imu) = deltasAhnfull(Nzz1,4)+i*deltasAhnfull(Nzz1,8);
  deltasT_Ahn_full(imu)   = deltasAhnfull(Nzz1,9)+i*deltasAhnfull(Nzz1,10);
end
%% This plot proves when phase=pi/4, the line mu=0 becomes the symmetry axis.
%% For mu'=-mu, Re(delta(mu'))=Im(delta(mu)), Im(delta(mu'))=Re(delta(mu)).
plot(mufull, real(deltasc_Ahn_full(:)), mufull, imag(deltasc_Ahn_full(:)))
plot(mufull, abs(deltasc_Ahn_full).^2)
%% just to check symmetry-------------------------------------------------------end


