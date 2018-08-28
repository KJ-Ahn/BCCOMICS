%%%% Script: Calculates (scale factor) ~ (growth factor) tables:
%%           azz~{Dplus_grow, Dplus_decay, Dminus_stream, 
%%                dDplus_grow_da, dDplus_decay_da, dDminus_stream_da}
%%
%% Author: Kyungjin Ahn
%%
%%
%% IMPORTANT: Accuracy of this script far surpasses what has been
%% described in Appendix A of A16. If interested in the details,
%% read comments in this script.
%%
%% It turns out that getting the asymptotes correctly is the key to
%% accuracy of the calculated growth factors. In the presence of
%% OmLambda0 and Omr0, the analytic form for growth factors is never
%% a good approximation unless a<<1e-3 (or z>>1000).
%%
%% "Fdecay" and "Fstream" equations try to pick up, unintendedly,
%% growing and compensated modes respectively, if integrated from
%% very high z(~10^7) such that asymptotes are the simplest (purely
%% radiation dominated). "Fgrow" just works fine starting from z=10^7.
%% So for "Fdecay" and "Fstream", one needs a different starting redshift.
%% This then requires one to integrate forward, and backward, and
%% stitch them afterwards.
%%
%% To find the asymptotes, I find it most optimal to have a 2nd-order
%% accurate expansion in polynomials, and integrating from where
%% epsr~3*epsL, where epsr=Omr0/Omega0/a and epsL=OmLambda0/Omega0*a^3.
%% At z=12, epsr~3*epsL, so z=12 becomes the starting point of integration.
%% Omega_matter = (1+epsr+epsL)^-1 ~ 1-(epsr+epsL)+(epsr+epsL)^2
%% H = H0*Om0^(-3/2)*(1+epsr+epsL)^(1/2) ~ expand to 2nd order.
%% F = 1+eps+eta, where eps and eta are 1st and 2nd order solutions,
%% respectively. After working tediously on evolution equations, one
%% can find eps and eta in terms of epsr and epsL.
%%
%% The reason why epsr~3*epsL is preferred is because evolution ODEs
%% are most sensitive to dF/da, and d(epsL)/da = 3* d(epsr)/da while
%% d(esr)/da = -epsr/a. So it is best to balance these two in dF/da.
%%
%% As a sanity check, one can do the same split-integration on "Fgrow"
%% and compare to one-swoop integration of "Fgrow". They match nicely.
%%%%

%% designated scale factors for a-D tables.
log10a      = linspace(log10(ai), log10(1), 100000)';
azz         = 10.^log10a;
log10az_min = log10a(1);
dlog10az    = log10a(2)-log10a(1);

%% ODE calculation option
options           = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% -- growing mode ------------ default ------------------------ begin
%% -- asymptote at a<<1e-3: Fgrow \propto a^-1, dFgrow/da = -Fgrow/a
%% -- With this small-a asymptote, integrating the 2nd order ODE 
%% -- for Fgrow correctly picks up the growing mode, so this works.
%% -- Starting Fgrow can be any value; normalization will follow.
abegin            = 1e-7;
Fgrow_init        = 1e3;
F0                = [Fgrow_init; -Fgrow_init/abegin];
options           = odeset('RelTol',1e-10,'AbsTol',1e-10);
[agode, Fgrowode] = ode45(@Fgrow, [abegin, 1], F0, options);
Fgrowode1         = Fgrowode(:,1);
Fgrowode2         = Fgrowode(:,2);
%% -- growing mode ------------ default ------------------------ end

%%%% -- growing mode ------------ optional ----------------------- begin
%%%% -- Calculating by backward and forward integration, if to make
%%%% -- comparison with the above "default" integration of Fgrow.
%%abegin              = 1/(1+12);
%%Fgrow_init          = 1;
%%epsr                = Omr0/Om0/abegin;
%%epsL                = OmLambda0/Om0*abegin^3;
%%dFda_init           = Fgrow_init/abegin*(-(2/3)*epsr-(6/11)*epsL+(96/187)*epsL^2+(64/99)*epsr*epsL)/(1+(2/3)*epsr-(2/11)*epsL+(16/187)*epsL^2+(32/99)*epsr*epsL);
%%F0                  = [Fgrow_init; dFda_init];
%%[agodeB, FgrowodeB] = ode45(@Fgrow, [abegin, 1/(1+1000)], F0, options);
%%FgrowodeB1          = FgrowodeB(:,1);
%%FgrowodeB2          = FgrowodeB(:,2);
%%
%%[agodeF, FgrowodeF] = ode45(@Fgrow, [abegin, 1], F0, options);
%%FgrowodeF1          = FgrowodeF(:,1);
%%FgrowodeF2          = FgrowodeF(:,2);
%%
%%%% -- stitch backward & forward integrations
%%agode               = [flipud(agodeB    );     agodeF(2:length(agodeF),1)];
%%Fgrowode1           = [flipud(FgrowodeB1); FgrowodeF1(2:length(agodeF),1)];
%%Fgrowode2           = [flipud(FgrowodeB2); FgrowodeF2(2:length(agodeF),1)];
%%%% -- growing mode ------------ optional ----------------------- end

%% -- decaying mode
%% -- Starting Fdecay can be any value; normalization will follow.
%% --
%% -- This asymptote for dFdecay/da (dFda_init) is correct to 2nd order.
abegin               = 1/(1+12); %% z=12
Fdecay_init          = 1;
epsr                 = Omr0/Om0/abegin;
epsL                 = OmLambda0/Om0*abegin^3;
dFda_init            = Fdecay_init/abegin*((9/14)*epsr-(25/28)*epsr^2+(3/2)*epsL-(3/4)*epsL^2+(9/14)*epsr*epsL)/(1-(9/14)*epsr+(25/56)*epsr^2+(1/2)*epsL-(1/8)*epsL^2+(9/28)*epsr*epsL);
F0                   = [Fdecay_init; dFda_init];
[adodeB, FdecayodeB] = ode45(@Fdecay, [abegin, 1/(1+1000)], F0, options);
FdecayodeB1          = FdecayodeB(:,1);
FdecayodeB2          = FdecayodeB(:,2);
%% -- Need to backward-integrate, but some old versions of ode45 do not allow this.
%% -- In this case, integration halts and the output array is just F0.
%% -- So stop the code in this case, and you need to upgrade the ode45 version.
if (length(adodeB) == 1)
  disp('Installed ode45 version not allowing backward integration; Halting.');
  disp('Please update ode package for gnu octave.');
  return;
end

[adodeF, FdecayodeF] = ode45(@Fdecay, [abegin, 1], F0, options);
FdecayodeF1 = FdecayodeF(:,1);
FdecayodeF2 = FdecayodeF(:,2);

%% -- stitch backward & forward integrations
adode      = [flipud(adodeB     );      adodeF(2:length(adodeF),1)];
Fdecayode1 = [flipud(FdecayodeB1); FdecayodeF1(2:length(adodeF),1)];
Fdecayode2 = [flipud(FdecayodeB2); FdecayodeF2(2:length(adodeF),1)];

%% -- streaming mode
%% -- Starting Fstream can be any value; normalization will follow.
%% --
%% -- This asymptote for dFstream/da (dFda_init) is correct to 2nd order.
abegin       = 1/(1+12); %% z=12
Fstream_init = 1;
epsr         = Omr0/Om0/abegin;
epsL         = OmLambda0/Om0*abegin^3;
dFda_init    = Fstream_init/abegin*((1/6)*epsr-(3/20)*epsr^2+(3/10)*epsL-(9/44)*epsL^2-(1/2)*epsr*epsL)/(1-(1/6)*epsr+(3/40)*epsr^2+(1/10)*epsL-(3/88)*epsL^2-(1/4)*epsr*epsL);
F0           = [Fstream_init; dFda_init];

[asodeB, FstreamodeB] = ode45(@Fstream, [abegin, 1/(1+1000)], F0, options);
FstreamodeB1 = FstreamodeB(:,1);
FstreamodeB2 = FstreamodeB(:,2);

[asodeF, FstreamodeF] = ode45(@Fstream, [abegin, 1], F0, options);
FstreamodeF1 = FstreamodeF(:,1);
FstreamodeF2 = FstreamodeF(:,2);

%% -- stitch backward & forward integrations
asode       = [flipud(asodeB      );       asodeF(2:length(asodeF),1)];
Fstreamode1 = [flipud(FstreamodeB1); FstreamodeF1(2:length(asodeF),1)];
Fstreamode2 = [flipud(FstreamodeB2); FstreamodeF2(2:length(asodeF),1)];


%% Now, multiply the analytic forms in completely matter-dominated epoch to
%% the 'F' factors to get the final growth factors 'D', then interpolate at
%% designated scale factors. See a few lines just above Eq. 26 in A16.
%% 
%% Best to interpolate in log space.
%%
%% Normalization to make F=1 @ z=1000 will follow.
log10Dplus_grow   = interp1(log10(agode), log10(Fgrowode1.*agode          ),log10a,'spline');
Dplus_grow        = 10.^log10Dplus_grow ;
dDplus_grow_da    = interp1(log10(agode),       Fgrowode2.*agode+Fgrowode1 ,log10a,'spline');

log10Dplus_decay  = interp1(log10(adode), log10(Fdecayode1.*adode.^-1.5                            ),log10a, 'spline');
Dplus_decay       = 10.^log10Dplus_decay;
dDplus_decay_da   = interp1(log10(adode),       Fdecayode2.*adode.^-1.5-1.5*Fdecayode1.*adode.^-2.5 ,log10a,'spline');

log10Dminus_stream = interp1(log10(asode), log10(Fstreamode1.*asode.^-0.5                             ),log10a,'spline');
Dminus_stream      = 10.^log10Dminus_stream;
dDminus_stream_da  = interp1(log10(asode),       Fstreamode2.*asode.^-0.5-0.5*Fstreamode1.*asode.^-1.5 ,log10a,'spline');


%% rescale to be 1 at zi(=1000). differentiated ones should use the same normalization coefficient.
dDplus_grow_da    = dDplus_grow_da   /Dplus_grow   (1);
dDplus_decay_da   = dDplus_decay_da  /Dplus_decay  (1);
dDminus_stream_da = dDminus_stream_da/Dminus_stream(1);

Dplus_grow        = Dplus_grow   /Dplus_grow   (1);
Dplus_decay       = Dplus_decay  /Dplus_decay  (1);
Dminus_stream     = Dminus_stream/Dminus_stream(1);

Fplus_grow_normed    = Dplus_grow    ./(azz/ai);
Fplus_decay_normed   = Dplus_decay   ./(azz/ai).^-1.5;
Fminus_stream_normed = Dminus_stream ./(azz/ai).^-0.5;

%% Output growth factors. This makes Fig 8. of A16. Now with increased accuracy,
%% F^d and D^d in Fig 8. of A16 has changed quantitatively.
gdata = [azz(1:10:length(azz)) Dplus_grow(1:10:length(azz)) Dplus_decay(1:10:length(azz)) Dminus_stream(1:10:length(azz)) Fplus_grow_normed(1:10:length(azz)) Fplus_decay_normed(1:10:length(azz)) Fminus_stream_normed(1:10:length(azz))];
fout = fopen([outputdir '/a_growth.dat'],'w');
fprintf(fout, '%e %e %e %e %e %e %e\n', gdata');
fclose(fout);
