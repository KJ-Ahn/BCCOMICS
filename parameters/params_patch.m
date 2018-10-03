%% Initial condition configuration inside a chosen patch
Lbox_p    = 1/h;         %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p   = 524;         %% # of cells & particles along one axis: make it an even number

ICdir     = '../ICs';  %% directory to place initial condition outputs

%% Nseed    denotes [setupdir '/subgaussseed' num2str(Nmode_p)  '.matbin'].
%% Noldseed denotes [setupdir '/subgaussseed' num2str(Noldseed) '.matbin'].
oldseedflag    = true;   %% if true, use old seed
%% Following two parameters needed only when oldseedflag = true
diroldseed     = setupdir %% directory where old seed exists
Noldseed       = 512     %% seed to use has Noldseed*Noldseed*(Noldseed/2+1) elements
%% When Noldseed = Ncell_p, the old seed is just the right choice.
%% When Noldseed > Ncell_p, part of the old seed will be used.
%% When Noldseed < Ncell_p, the old seed will be used and missing high-k seeds 
%% will be generated and attached.

%% name of old seed file = ['subgaussseed' num2str(Noldseed) '.matbin']

recordseedflag = true;   %% if true, record the used seed with right dimension
%% If following flag is true, Eulerian velocity field is interpolated at 
%% displaced particle positions for particle velocity, surpassing accuracy of 
%% 1st-order Lagrangian perturbation to some extent. 
%% If following flag is true memory usage will increase quite significantly.
particlevelocity_accuracyflag = true;  



