%% Calculate the initial("_i") and final("_f") cosmological parameters
%% for a given over(under)dense patch. The local cosmological values
%% are denoted by "_loc".

clear;
more off;

%% Read in constants in cgs unit and conversion factors.
Consts_Conversions;  %%==== script ==================
%% Read in parameters
run('params.m');  %%==== script ==================

%% Read in cosmology
run(Cosmology);  %%==== script ==================

fb = ombh2/(ombh2+omch2); %% baryon/matter fraction
fc = omch2/(ombh2+omch2); %% CDM/matter fraction

%% Read in parameters for initial condition
run('params_patch.m');  %%==== script ==================

global H_i H_l_i Om_l_i Omr_l_i OmLambda_l_i OmK_l_i aloci;

Lbox_p_inMpch = Lbox_p*h;  %% enzo uses 'ComovingBoxSize' in units of Mpc/h

mmw = 1.2195 %% mean molecular weight with X=0.76, Y=0.24, neutral.
rhocrit0 = 3*(100*h*1e5/Mpc)^2/(8*pi*G) %% in g/cm^3

%% read in the redshift. zzbegin=1000, and zzend is for the
%% redshift of the initial conditions.
zz    = load('zz.dat');
zzbegin = zz(1);
zzend   = zz(2);
azbegin = 1/(1+zzbegin);
azend   = 1/(1+zzend);

iccdat = load('icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat');
icc = iccdat(1:3);
Dc  = iccdat(4);
Db  = iccdat(5);
Thc = iccdat(6);

fin = fopen('stats_zend.dat','r');
dum = fgetl(fin);
dum = fgetl(fin);
datstat = fscanf(fin,'%e %e %e %e %e %e %e %e %e')';
fclose(fin);

flagmean = (abs(Dc)<1e-3); %% flag for mean-density cases

H_loc_i         = zeros(Ncc,1);
rhocrit_loc_i   = zeros(Ncc,1);
rhocrit_ratio_i = zeros(Ncc,1);
Om_loc_i        = zeros(Ncc,1);
Omr_loc_i       = zeros(Ncc,1);
OmLambda_loc_i  = zeros(Ncc,1);
OmK_loc_i       = zeros(Ncc,1);

%% For the same-cosmic-time Enzo outputs, first list wanted global redshifts
zglobal_enzo = load('zglobal.dat');
zglobal_enzo = sort(zglobal_enzo, 'descend'); %% sort in descending order
%zglobal_enzo = [linspace(200,120,5)';  linspace(100,40,7)'; linspace(30,22,5)'; linspace(20,11,10)'; linspace(10,3,15)']
aglobal_enzo = 1./(1+zglobal_enzo);
Nz_enzo      = length(zglobal_enzo);

thefactor = sqrt(Om0/azend^3 + Omr0/azend^4 + OmLambda0); % global one

%% lock fb_l/fc_l ratio locked for the patch. In practice very close to fb/fc.
fc_l         = (1+Dc)*fc / ((1+Dc)*fc + (1+Db)*fb); %% local CDM fraction
fb_l         = 1 - fc_l; %% local CDM fraction
Ddot_over_D1 = -Thc/(1+Dc);  %% Myr^-1, using dD/dt=-Th relation, and follow CDM only.

H_i          = H0*thefactor; %% initial Hubble constant (Myr^-1) for global, flat universe
H_loc_i = H_i - (1/3)*Ddot_over_D1; %% As in Goldberg & Vogeley (2004, eq. 3)

rhocrit_i       = 3*(H_i    *s_inMyr)^2 / (8*pi*G); %% g/cm^3
rhocrit_loc_i   = 3*(H_loc_i*s_inMyr)^2 / (8*pi*G); %% g/cm^3
rhocrit_ratio_i = rhocrit_loc_i/rhocrit_i; 

%% initial global Omega's.
Om_i       = (Om0 /azend^3) / thefactor^2;
Omr_i      = (Omr0/azend^4) / thefactor^2;
OmLambda_i = (OmLambda0)    / thefactor^2;

%% initial local Omega's.
Om_loc_i       = Om_i*(1+Dc_tab)/ rhocrit_ratio_i;
Omr_loc_i      = Omr_i          / rhocrit_ratio_i;
OmLambda_loc_i = OmLambda_i     / rhocrit_ratio_i;
OmK_loc_i      = 1 - (Om_loc_i + Omr_loc_i + OmLambda_loc_i);

%% record useful values
fout        = fopen('global_and_local_quantities.dat','w');
datglobal_i = [zzend rhocrit_i H_i Om_i Omr_i OmLambda_i];
fprintf(fout, '## zred rhocrit(g/cm^3) H(Myr^-1) Om Omr OmLambda\n');
fprintf(fout, '%i %e %e %e %e %e\n', datglobal_i');
fprintf(fout, '\n');

datloc_i = [icc_tab rhocrit_loc_i H_loc_i Om_loc_i Omr_loc_i OmLambda_loc_i OmK_loc_i];
fprintf(fout, '## icc1 icc2 icc3 rhocrit(g/cm^3) H(Myr^-1)  Om   Omr   OmLambda  OmK\n');
fprintf(fout, '%4i %4i %4i %e %e %e %e %e %e\n', datloc_i');
fclose(fout);


%%%%%%%%% Now do integration to obtain (time) ~ (scale factor) table,
%%%%%%%%% both for the global case and for each patch.

%%%% for global case
a_i       = azend   %% azend corresponds to "_i" as explained above.
TimeUnits = 2.519445e17/sqrt(Om0) /h * a_i^(3/2)  %% Enzo time unit (in seconds), which is related to a_i. See CosmologyGetUnits.C from enzo src. The number seems to be using a slightly different Myr calculation.

%% For simplicity, set the initial time = 0.
%% Whenever required, one can add the actual cosmic time_i*H_i to the timetable.

%% assign values to globals: mean-density patch just follows global LCDM evolution
H_l_i        = H_i;
Om_l_i       = Om_i;
Omr_l_i      = Omr_i;
OmLambda_l_i = OmLambda_i;
OmK_l_i      = 1 - (Om_l_i + Omr_l_i + OmLambda_l_i);
aloci        = a_i;

%%%% time ~ (scale factor) table for global case.
%% take small enough value for radiation domination
tiHi                 = 0.000001; 
tfHi                 = 1000;
options              = odeset('RelTol',1e-6,'AbsTol',1e-9);
%% initial a value, assuming radiation domination, is given analytically.
[tHiode, aglobalode] = ode45(@fdadt, [tiHi, tfHi], sqrt(2*tiHi*sqrt(Omr_i))*a_i, options);

tHi_a_0p25     = interp1(aglobalode, tHiode, 0.25, 'spline');
%% tHi table corresponding to aglobal table.
tHiglobal_enzo = interp1(aglobalode, tHiode, aglobal_enzo, 'spline');

loglog(tHiode, aglobalode, tHiode, 6.4e-3*tHiode.^(2/3))
axis([1e-3 1e4 5e-5 1])

dattemp = [tHiode aglobalode];
fout=fopen('tHi_a_fine.dat','w');  %% columns: t*Hi, a. (Hi is again H at given initial time)
fprintf(fout,'%e %e\n', dattemp');
fclose(fout);

dattemp = [tHiglobal_enzo aglobal_enzo];
fout=fopen('tHi_a.dat','w');  %% columns: t*Hi, a. (Hi is again H at given initial time)
fprintf(fout,'%e %e\n', dattemp');
fclose(fout);

%%%% time ~ (scale factor) table for local case.

for icc=1:Ncc
    H_l_i        = H_loc_i(icc);
    Om_l_i       = Om_loc_i(icc);
    Omr_l_i      = Omr_loc_i(icc);
    OmLambda_l_i = OmLambda_loc_i(icc);
    OmK_l_i      = 1 - (Om_l_i + Omr_l_i + OmLambda_l_i);
    aloci        = a_i; %% To make Enzo use the same code units, a_i should be universal

    %% Check if turnaround occurs (not completed yet!!!)
    aloc_test    = linspace(a_i,1,1e5);
    insidesquare = Om_l_i./(aloc_test/aloci) + Omr_l_i./(aloc_test/aloci).^2 + OmLambda_l_i*(aloc_test/aloci).^2 + OmK_l_i;
    if (any(insidesquare<0))
        [xval,fval]=fzero(@(x) Om_l_i/x+Omr_l_i/x^2+OmLambda_l_i*x^2+OmK_l_i, [1,1/aloci]) 
        a_loc_turnaround = xval*aloci
        %% Do something clever to change tfHi
    end
    
    tfHi = tHi_a_0p25; %% Choose universal end time, which is "0" or "present" for each icc.
    [tHiode_l, aglobalode_l] = ode45(@fdadt, [tiHi, tfHi], sqrt(2*tiHi*H_l_i/H_i*sqrt(Omr_l_i))*aloci, options);
    aloc_enzo(:,icc) = interp1(tHiode_l, aglobalode_l, tHiglobal_enzo, 'spline');
    alocf(icc,1) = aglobalode_l(length(aglobalode_l));

end

clear dattemp;
dattemp = [tHiglobal_enzo aloc_enzo];
fout=fopen('tHi_alocal.dat','w');  %% columns: t*Hi, alocal(icc). (Hi is again H at given initial time)
fprintf(fout,'%e %e %e %e %e %e %e %e %e %e %e %e %e\n', dattemp'); %% should have Ncc+1 columns.
fclose(fout);

%%%% for each alocf(icc,i), calculate cosmological parameters at that time.
for icc=1:Ncc
    aratio             = alocf(icc)/aloci;
    denominator(icc,1) = Om_loc_i(icc)*aratio^(-3) + Omr_loc_i(icc)*aratio^(-4) + OmLambda_loc_i(icc) + OmK_loc_i(icc)*aratio^(-2);
    Om0_l      (icc,1) = Om_loc_i(icc)      *aratio^(-3) / denominator(icc);
    Omr0_l     (icc,1) = Omr_loc_i(icc)     *aratio^(-4) / denominator(icc);
    OmLambda0_l(icc,1) = OmLambda_loc_i(icc)             / denominator(icc);
    OmK0_l     (icc,1) = OmK_loc_i(icc)     *aratio^(-2) / denominator(icc);
    h0_l       (icc,1) = H_loc_i(icc)/(100*km_inMpc/s_inMyr) * sqrt(denominator(icc));  %% unitless

    %% Enzo uses "0" values to denote when the scale factor = 1, and
    %% redshift = 0. In order to satisfy this convention AND make the time
    %% of alocf(icc) correspond to "0", from the relation 
    %% aloc/alocf = (1+zlocf)/(1+zloc) = 1/(1+zloc_new),
    %% 1+zloc_new = (1+zloc)/(1+zlocf), or zloc_new = (1+zloc)/(1+zlocf) -1.
    %% This also means                     zloc_new = alocf/aloc -1.
    %% Similarly, aloc/alocf = aloc_new.

    aloc_new_enzo(:,icc) = aloc_enzo(:,icc)/alocf(icc);
    zloc_new_enzo(:,icc) = 1./aloc_new_enzo(:,icc) - 1;
end

%%%% Write an enzo parameter file patch, to be included in there.
foutratio=fopen('enzounit_ratio.dat','w');
for icc=1:Ncc
    strbox = ['Lbox-' num2str(Lbox_p_inMpch) 'Mpch_'];
    stric = num2str(icc_tab(icc,1));
    strjc = num2str(icc_tab(icc,2));
    strkc = num2str(icc_tab(icc,3));
    strout = [strbox 'ic' stric '_jc' strjc '_kc' strkc '_enzoparm.dat'];
    fout = fopen(strout, 'w');
    fprintf(fout, 'CosmologySimulationOmegaBaryonNow        = %f\n', Om0_l(icc)*fb_l(icc) );
    fprintf(fout, 'CosmologySimulationOmegaCDMNow           = %f\n', Om0_l(icc)*fc_l(icc) );
    fprintf(fout, '\n');
    fprintf(fout, 'CosmologyOmegaMatterNow    = %f\n', Om0_l(icc)       );
    fprintf(fout, 'CosmologyOmegaLambdaNow    = %f\n', OmLambda0_l(icc) );
    fprintf(fout, 'CosmologyOmegaRadiationNow = %f\n', Omr0_l(icc)      );
    fprintf(fout, 'CosmologyHubbleConstantNow = %f\n', h0_l(icc)        );
    %% comoving box size is the proper size at "0", so it is (proper size at a_i)*(expansion ratio)
    %%                                                     = L * a_i * (alocf/aloci) = L * alocf
    %% But also, Lbox_p_inMpch uses h, not h0_l, so need to rescale with (h0_l/h)
    fprintf(fout, 'CosmologyComovingBoxSize   = %f  ', Lbox_p_inMpch*(h0_l(icc)/h)*alocf(icc)); 
    fprintf(fout, ' // Mpc/h\n' );  
    fprintf(fout, 'CosmologyInitialRedshift   = %f\n', zloc_new_enzo(1, icc)        );
    fprintf(fout, 'CosmologyFinalRedshift     = %f\n', zloc_new_enzo(Nz_enzo,icc)   );
    fprintf(fout, '\n');
    for iz_enzo = 1:Nz_enzo
        fprintf(fout, 'CosmologyOutputRedshift[%i] = %f\n', iz_enzo-1, zloc_new_enzo(iz_enzo,icc) );
    end
    fclose(fout);
    %% Print out the ratios of enzo units (in CosmologyGetUnits.C).
    %% ratio = global/local
	lengunit_ratio = Lbox_p_inMpch/(h*(1+zzend)) /(Lbox_p_inMpch*(h0_l(icc)/h)*alocf(icc)/(h0_l(icc)*(1+zloc_new_enzo(1,icc))));
    densunit_ratio = Om0*h^2*(1+zzend)^3 /(Om0_l(icc)*h0_l(icc)^2*(1+zloc_new_enzo(1,icc))^3);
    timeunit_ratio = sqrt(Om0_l(icc))*h0_l(icc)*(1+zloc_new_enzo(1,icc))^1.5 /(sqrt(Om0)*h*(1+zzend)^1.5);
    velounit_ratio = Lbox_p_inMpch*sqrt(Om0)*sqrt(1+zzend) /(Lbox_p_inMpch*(h0_l(icc)/h)*alocf(icc)*sqrt(Om0_l(icc))*sqrt(1+zloc_new_enzo(1,icc)));
    tempunit_ratio = velounit_ratio^2;
    datcc = [icc flagmean(icc) icc_tab(icc,:) lengunit_ratio densunit_ratio timeunit_ratio velounit_ratio tempunit_ratio];
    fprintf(foutratio, '%2i %2i  %4i %4i %4i  %e %e %e %e %e\n', datcc');
end
fclose(foutratio);











