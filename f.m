%% Growth rate equations for a given mode (in k-space) parameterized by
%% x(1)~x(10) inside a given patch (Eq. 11 of A16, Eq. 2 of AS18), but
%% in /da instead of /dt, where a is the scale factor.
%%
%% x(1)~x(10): real and imaginary values of delta_c, etc. See "What
%%             are {x}" below for details. 
%%
%% Values of patches-
%%   Delta, Theta, rV = V_cb = -Vbc: Eq. 7 of A16
%%   Delta_T:                        Eq. 9 of A16
%%


function dxda = f(a,x)
  global mH kb MpcMyr_2_kms H0 Om0 Omr0 OmLambda0 TCMB0 ai aa1 aa2 fb fc Thc_i Thb_i rV_i;
  global Dplus_grow Dplus_decay Dminus_stream dDplus_grow_da dDplus_decay_da azz log10az_min dlog10az;
  global Deltagro_kk Deltadec_kk Deltacom_kk Deltastr_kk;
  global signDT alpha coeff_Delta_T beta gamma;
  global log10az_min dlog10az
  global zrecf xerecf dzrecf zrecf1 tgamma;
  global ksample costh;

  dxda    = zeros(10,1);

  %% What are {x}:
  %% x(1)=Re(delta_c); x(2)=Re(theta_c); x(3)=Re(delta_b); x(4)=Re(theta_b)
  %% x(5)=Im(delta_c); x(6)=Im(theta_c); x(7)=Im(delta_b); x(8)=Im(theta_b)
  %% x(9)=Re(delta_T); x(10)=Im(delta_T)
  %% dx/da = (dx/dt)/(a*H)
  Hz      = H0*sqrt(Om0*a^-3+Omr0*a^-4 + OmLambda0);
  Omz     = Om0*a^-3 /(Om0*a^-3 + Omr0*a^-4 + OmLambda0);
  aH      = a*Hz;
  Tbz     = TCMB0/a /(1+a/aa1/(1+(aa2/a)^1.5)); %% fit by TH10.
  Tgamma  = TCMB0/a;
  %% 1e-10 required to convert cm/s to km/s, and get sound speed 
  %% in Mpc/Myr unit.
  por     = kb*Tbz/(1.22*mH) * 1e-10 / MpcMyr_2_kms^2; %% p over rho

  %% Interpolate growth factors without calling the whole array
  %% Works only when the base table is in uniform in log(a).
  ind1 = floor((log10(a)-log10az_min)/dlog10az)+1;
  ind2 = ind1+1;
  az1  = azz(ind1);
  az2  = azz(ind2);
  %% linear interpolation of growth factors at a given a.
  Dpg = (Dplus_grow(ind2)   -Dplus_grow(ind1)   )/(az2-az1)*(a-az1)+Dplus_grow(ind1); %% D^{g}(a) in A16, Eq. 7.
  Dpd = (Dplus_decay(ind2)  -Dplus_decay(ind1)  )/(az2-az1)*(a-az1)+Dplus_decay(ind1); %% D^{d}(a) in A16, Eq. 7.
  Dms = (Dminus_stream(ind2)-Dminus_stream(ind1))/(az2-az1)*(a-az1)+Dminus_stream(ind1); %% D^{s}(a) in A16, Eq. 7.
  dDpg_da = (dDplus_grow_da(ind2)-dDplus_grow_da(ind1))/(az2-az1)*(a-az1)+dDplus_grow_da(ind1); %% dD^{g}(a)/da in A16, Eq. 7.
  dDpd_da = (dDplus_decay_da(ind2)-dDplus_decay_da(ind1))/(az2-az1)*(a-az1)+dDplus_decay_da(ind1); %% dD^{d}(a)/da in A16, Eq. 7.

  %% now for large-scale Delta_c, Delta_b, Theta_c, Theta_b
  Delta_plus_g    = Deltagro_kk * Dpg;  %% Deltagro_kk = Delta_{gro} in A16, Eq. 7.
  Delta_plus_d    = Deltadec_kk * Dpd;  %% Deltadec_kk = Delta_{dec} in A16, Eq. 7.
  Delta_minus_str = Deltastr_kk * Dms;  %% Deltastr_kk = Delta_{str} in A16, Eq. 7.

  Theta_plus_g    = -aH* Deltagro_kk * dDpg_da;
  Theta_plus_d    = -aH* Deltadec_kk * dDpd_da;
  
  %% A16, Eq. 7.
  Delta_c = (Delta_plus_g + Delta_plus_d) +fb*(Deltacom_kk + Delta_minus_str);
  Delta_b = (Delta_plus_g + Delta_plus_d) -fc*(Deltacom_kk + Delta_minus_str);
  Theta_c = (Theta_plus_g + Theta_plus_d) +fb*(Thc_i-Thb_i)*(a/ai)^-2;
  Theta_b = (Theta_plus_g + Theta_plus_d) -fc*(Thc_i-Thb_i)*(a/ai)^-2;
  rV      = rV_i * (a/ai)^-1;

  %% Fitting formula for evolution of background temperature fluctuation
  %% when the Eulerian patch size is 4Mpc. Don't know for other size patches... 
  %% beta=2.8, gamma=0.33
  %% 0.2850509 = (log10(0.1)+2.8)^0.33 - (log10(0.01)+2.8)^0.33
  %%  alpha   = log10(DT_B/DT_A)/0.2850509: defined outside
  if (a<0.0018) %% just flatten out small Delta_T at z>~500.
    %% 0.38462 = (log10(0.0018)+2.8)^0.33
    %%    Delta_T = 10^(alpha*0.38462            ) * 10^(-Y);
    Delta_T = signDT * 10^(alpha*0.38462              ) * coeff_Delta_T;
  else
    %%    Delta_T = 10^(alpha*(log10(a)+2.8)^0.33) * 10^(-Y);
    Delta_T = signDT * 10^(alpha*(log10(a)+beta)^gamma) * coeff_Delta_T;
  end

  %% Linear interpolation of xe(z) from a table from recfast.
  zz   = 1/a -1;
  ind1 = floor((zrecf1-zz)/dzrecf)+1;
  ind2 = ind1+1;
  zz1  = zrecf(ind1);
  zz2  = zrecf(ind2);
  xe1  = xerecf(ind1);
  xe2  = xerecf(ind2);

  xez  = (xe2-xe1)/(zz2-zz1)*(zz-zz1)+xe1;

  %% take V_c=0 frame.
  %% A16, Eq. 11, but with real and complex numbers for small-scale k-space 
  %% fluctuation variables, such as delta_c.
  %% rV = -V_bc = V_cb = V_c - V_b. 

  dxda(1) = (                           -(1+Delta_c)*x(2)        -Theta_c*x(1)) /aH;
  dxda(2) = (                           -1.5*Hz^2*Omz*(fc*x(1)+fb*x(3)) -2*Hz*x(2)) /aH;
  dxda(3) = (-1/a*rV*ksample*costh*x(7) -(1+Delta_b)*x(4)        -Theta_b*x(3)) /aH;
  dxda(4) = (-1/a*rV*ksample*costh*x(8) -1.5*Hz^2*Omz*(fc*x(1)+fb*x(3)) -2*Hz*x(4) +por*ksample^2/a^2*((1+Delta_T)*x(3)+(1+Delta_b)*x(9) )) /aH;
  dxda(5) = (                           -(1+Delta_c)*x(6)        -Theta_c*x(5)) /aH;
  dxda(6) = (                           -1.5*Hz^2*Omz*(fc*x(5)+fb*x(7)) -2*Hz*x(6)) /aH;
  dxda(7) = ( 1/a*rV*ksample*costh*x(3) -(1+Delta_b)*x(8)        -Theta_b*x(7)) /aH;
  dxda(8) = ( 1/a*rV*ksample*costh*x(4) -1.5*Hz^2*Omz*(fc*x(5)+fb*x(7)) -2*Hz*x(8) +por*ksample^2/a^2*((1+Delta_T)*x(7)+(1+Delta_b)*x(10))) /aH;

  dxda(9)  = 2/3*(dxda(3)*(1+Delta_T-Delta_b) -Theta_b*(x(9) -x(3))/aH) - xez/tgamma/a^4*Tgamma/Tbz*x(9) /aH ; 
  dxda(10) = 2/3*(dxda(7)*(1+Delta_T-Delta_b) -Theta_b*(x(10)-x(7))/aH) - xez/tgamma/a^4*Tgamma/Tbz*x(10)/aH ; 
end
