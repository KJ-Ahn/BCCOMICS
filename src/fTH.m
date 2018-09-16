%% Growth rate equations by TH for a given mode (in k-space) parameterized by
%% x(1)~x(10) inside a given patch (Eq. 11 of A16, Eq. 2 of AS18), but
%% in /da instead of /dt, where a is the scale factor, and rV_i is the only
%% patch quantity in the equation.

function dxda = fTH(a,x)
  global mH kb MpcMyr_2_kms H0 Om0 Omr0 OmLambda0 TCMB0 ai aa1 aa2 fb fc Thc_i Thb_i rV_i;
  global zrecf xerecf dzrecf zrecf1 tgamma;
  global ksample costh;
  
  dxda = zeros(10,1);
  %% What are {x}:
  %% x(1)=Re(delta_c); x(2)=Re(theta_c); x(3)=Re(delta_b); x(4)=Re(theta_b)
  %% x(5)=Im(delta_c); x(6)=Im(theta_c); x(7)=Im(delta_b); x(8)=Im(theta_b)
  %% x(9)=Re(delta_T); x(10)=Im(delta_T)
  %% dx/da = (dx/dt)/(a*H)
  Hz      = H0*sqrt(Om0*a^-3+Omr0*a^-4 + OmLambda0);
  Omz     = Om0*a^-3 /(Om0*a^-3 + Omr0*a^-4 + OmLambda0);
  aH      = a*Hz;
  Tbz     = TCMB0/a /(1+a/aa1/(1+(aa2/a)^1.5));
  Tgamma  = TCMB0/a;
  %% 1e-10 required to convert cm/s to km/s, and get sound speed 
  %% in Mpc/Myr unit.
  por     = kb*Tbz/(1.22*mH) * 1e-10 / MpcMyr_2_kms^2; %% p over rho
  rV      = rV_i * (a/ai)^-1;

  %% interpolation of xe(z) from a table, done efficiently without calling
  %% interp1 with the whole array of zrecf & xerecf.
  zz       = 1/a -1;
  ind1    = floor((zrecf1-zz)/dzrecf)+1;
  ind2    = ind1+1;
  zz1     = zrecf(ind1);
  zz2     = zrecf(ind2);
  xe1     = xerecf(ind1);
  xe2     = xerecf(ind2);
  %% Use straightforward expression rather than calling interp1 (time-consuming)
  %% xez  = interp1([zz1;zz2],[xe1;xe2],zz,'linear','extrap');
  xez     = (xe2-xe1)/(zz2-zz1)*(zz-zz1)+xe1;

  %% take V_c=0 frame.
  dxda(1) = (                           -            x(2)                     ) /aH;
  dxda(2) = (                           -1.5*Hz^2*Omz*(fc*x(1)+fb*x(3)) -2*Hz*x(2)) /aH;
  dxda(3) = (-1/a*rV*ksample*costh*x(7) -            x(4)                     ) /aH;
  dxda(4) = (-1/a*rV*ksample*costh*x(8) -1.5*Hz^2*Omz*(fc*x(1)+fb*x(3)) -2*Hz*x(4) +por*ksample^2/a^2*(x(3)+x(9))) /aH;
  dxda(5) = (                           -            x(6)                     ) /aH;
  dxda(6) = (                           -1.5*Hz^2*Omz*(fc*x(5)+fb*x(7)) -2*Hz*x(6)) /aH;
  dxda(7) = ( 1/a*rV*ksample*costh*x(3) -            x(8)                     ) /aH;
  dxda(8) = ( 1/a*rV*ksample*costh*x(4) -1.5*Hz^2*Omz*(fc*x(5)+fb*x(7)) -2*Hz*x(8) +por*ksample^2/a^2*(x(7)+x(10))) /aH;
  dxda(9)  = 2/3*dxda(3) - xez/tgamma/a^4*Tgamma/Tbz*x(9) /aH ; %% stripped down eqtn for dxda(9), real part
  dxda(10) = 2/3*dxda(7) - xez/tgamma/a^4*Tgamma/Tbz*x(10)/aH ; %% stripped down eqtn for dxda(9), imaginary part
end

