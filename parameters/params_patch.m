%% Initial condition configuration inside a chosen patch
Lbox_p    = 1/h;         %% box size in Mpc unit; Should be <= Lbox/Ncell
Ncell_p   = 508;         %% # of cells & particles along one axis: make it an even number

ICdir     = '../ICs';  %% directory to place initial condition outputs

%% Nseed   denotes [setupdir '/subgaussseed' num2str(Nmode_p) '.matbin'].
%% 512seed denotes [setupdir '/subgaussseed512.matbin'].
oldseedflag    = true;   %% if true, use old seed
%% Following two parameters needed only when oldseedflag = true
diroldseed     = setupdir %% directory where old seed exists
Noldseed       = 512     %% seed with Noldseed*Noldseed*(Noldseed/2+1) elements
%% name of old seed file = ['subgaussseed' num2str(Noldseed) '.matbin']

recordseedflag = true;   %% if true, record newly generated seed
%% record512companionflag will be automatically set to false when oldseedflag = true.
%% The following, if true, will overwrite existing 512seed.
record512companionflag = false; %% if true, generate companion 512seed and record it. 


