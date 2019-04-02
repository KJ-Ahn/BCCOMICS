%% Initial condition configuration inside a chosen patch
Lbox_p    = 1/h;         %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p   = 64;         %% # of cells & particles along one axis: make it an even number

ICdir     = '../ICs';  %% directory to place initial condition outputs

%% Nseed    denotes [setupdir '/subgaussseed' num2str(Nmode_p)  '.matbin'].
%% Noldseed denotes [setupdir '/subgaussseed' num2str(Noldseed) '.matbin'].
oldseedflag    = false;   %% if true, use old seed; if false, generate new seed
%% Following two parameters needed only when oldseedflag = true
diroldseed     = '/home/kjahn/BCCOMICS/ICs/1.00Mpch_64_ic40_jc84_kc79'; %% directory where old seed exists
%%diroldseed     = setupdir; %% directory where old seed exists
Noldseed       = 64;     %% seed to use has Noldseed*Noldseed*(Noldseed/2+1) elements
%% When Noldseed = Ncell_p, the old seed is just the right choice.
%% When Noldseed > Ncell_p, part of the old seed will be used.
%% When Noldseed < Ncell_p, the old seed will be used and missing high-k seeds 
%% will be generated and attached.

%% Name of old seed file = [diroldseed '/subgaussseed' num2str(Noldseed) '.matbin']

baryonparticleflag = true; %% if true, record baryon particle data (position & velocity).
recordseedflag = true;   %% if true, record the used seed with right dimension

%% If following flag is true, Eulerian velocity field is interpolated at 
%% displaced particle positions for particle velocity, surpassing accuracy of 
%% 1st-order Lagrangian perturbation theory (1LPT) to some extent. 1LPT just
%% uses initial Eulerian density to get velocity field, so 1LPT makes particle
%% (displaced from cell center) velocity equal to cell velocity, which is 
%% obviously poor in accuracy.
%% Same interpolation is done also for thermal energy of baryon particles.
%% If following flag is set to true, memory usage will increase though.
particlevelocity_accuracyflag = true;


%% Choose which simulation code & which format (may overlap)
enzo_bin_flag   = true;
enzo_HDF5_flag  = true;
gadget_bin_flag = false;  %% Not implemented yet
