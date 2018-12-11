%% This script generates initial condition data for enzo and dumps them.
%% For ICs for other simulation codes, this is the script to start from.
%%
%% Initial condition generation using following steps, for enzo.
%% (1) Get k and mu dependent fluctuation using transfer function.
%% (2) Apply random seed and normalize (rand_real_norm).
%% (3) FFT & record.

enzo_HDF5_flag = (enzo_HDF5_flag && matlabflag);  %% HDF5 output only possible in Matlab
%% Output binary when running Octave & the user wrongfully intends HDF5 output
if (~enzo_bin_flag && ~matlabflag)  
  disp('======= Octave cannot write HDF5 files for enzo. ');
  disp('======= Will instead write plain binary files. ');
  disp('======= You need to do "python convert2enzo.py" after copying');
  disp('======= convert2enzo.py from BCCOMICS/src_converter/ to the');
  disp(['======= directory' ICsubdir]);
  enzo_bin_flag = true;
end

%% =========== CDM density and position ======================== begin

%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
%% with array size Nmode_p, Nmode_p, Nmode_p.
%% Matlab allows extrapolation only for 'spline' method.
%% (costh_k_V and ksq_p have same array dimension, so interp2 
%%  interprets these as scattered data points: see interp2 instruction)
%% This is linear logarithmic interpolation along k, so the monopole
%% term (k=0) may obtain inf or nan due to 0.5*log(ksq_p). 
%% We will cure this by nullifying monopole anyway down below (**).
disp('----- Interpolating transfer function -----');
dc  = interp2(muext,log(ksampletab), deltasc,  costh_k_V,0.5*log(ksq_p),interp2opt);  %% dc still k-space values here.

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
dc = rand_real_norm(dc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);


%% CDM displacement vector, related to CDM density at 1st order.
%% No need for above normalization because this is
%% derived after above normalization on dc.
%% ------------- cpos1 ----------------------
disp('----- Calculating CDM position x -----');
Psi1                 = i*k1_3D_p./ksq_p.*dc;
Psi1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi1                 = real(ifftn(ifftshift(Psi1)));

xCDM_plane    =   Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
xCDM_ex_plane = 5*Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in CDM position

%% Normalized position of particles in domain [0,1), cell-centered way. (enzo)
%% If unperturbed(Psi=0), it should run [0.5, 1.5, ...., Nmode_p-0.5]/Nmode_p,
%% For enzo, wrapping needed if perturbed potition is out of the domain [0, 1).
Psi1 = mod((Psi1 + (k1_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);

if enzo_bin_flag
  fout = fopen([ICsubdir '/cpos1'], 'w');
  fwrite(fout, Psi1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Prepare to write Particle Positions in hdf5
  datname     = 'ParticlePositions';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  ND1 = Ncell_p;
  ND3 = Ncell_p^3;
  h5create(foutname,datasetname,[ND3 3]);
  h5write(foutname,datasetname,Psi1(:), [1 1], [ND3 1]);
end

%% Prepare for interpn, let positions run from 1:Nmode_p for unperturbed particles,
%% to get more-accurate-than-1LPT velocity when particlevelocity_accuracyflag=true.
%% Using mod function, the actual positions will run from 1 to Nmode_p+0.9999999...
%% Interpolation basis will have domain 1:Nmode_p+1 for safe interpolation (see **).
if particlevelocity_accuracyflag  
  Psi1 = mod(Psi1*Nmode_p - 0.5, Nmode_p)+1;
else
  clear Psi1  %% save memory
end

%% ------------- cpos2 ----------------------
disp('----- Calculating CDM position y -----');
Psi2                 = i*k2_3D_p./ksq_p.*dc;
Psi2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi2                 = real(ifftn(ifftshift(Psi2)));

yCDM_plane    =   Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
yCDM_ex_plane = 5*Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in CDM position

Psi2 = mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);

if enzo_bin_flag
  fout = fopen([ICsubdir '/cpos2'], 'w');
  fwrite(fout, Psi2, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,Psi2(:), [1 2], [ND3 1]);
end

if particlevelocity_accuracyflag
  Psi2 = mod(Psi2*Nmode_p - 0.5, Nmode_p)+1;
else
  clear Psi2  %% save memory
end

%% ------------- cpos3 ----------------------
disp('----- Calculating CDM position z -----');
Psi3                 = i*k3_3D_p./ksq_p.*dc;
Psi3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi3                 = real(ifftn(ifftshift(Psi3)));

zCDM_plane    =   Psi3(:,:,1) + k3_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
zCDM_ex_plane = 5*Psi3(:,:,1) + k3_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in CDM position

Psi3 = mod((Psi3 + (k3_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);

if enzo_bin_flag
  fout = fopen([ICsubdir '/cpos3'], 'w');
  fwrite(fout, Psi3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,Psi3(:), [1 3], [ND3 1]);
  %% Finalize writing Particle Positions in hdf5
  topgriddims = -99999*ones(1,3);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(3));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(1));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(ND3));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims-1));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end

if particlevelocity_accuracyflag
  Psi3 = mod(Psi3*Nmode_p - 0.5, Nmode_p)+1;
else
  clear Psi3  %% save memory
end

%% ------------- density ---------------------
dc = real(ifftn(ifftshift(dc)));  %% just for debugging

Zc = reshape(dc(:,:,1),Nmode_p,Nmode_p); %% for figure
clear dc  %% save memory
%% =========== CDM density and position ======================== end



%% =========== CDM velocity ==================================== begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
disp('----- Interpolating transfer function -----');
Thc = interp2(muext,log(ksampletab), deltasThc, costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
Thc = rand_real_norm(Thc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vc1 ----------------------
disp('----- Calculating CDM velocity x -----');
vc1(:,:,:)          = -i*af*k1_3D_p./ksq_p.*Thc(:,:,:);
vc1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc1                 = real(ifftn(ifftshift(vc1)));

Vc1 = reshape(vc1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity x more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc1 = padarray(vc1, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc1 = interpn(vc1, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
%% into enzo velocity unit
vc1 = vc1 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vc1'], 'w');
  fwrite(fout, vc1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Prepare to write Particle Velocities in hdf5
  datname     = 'ParticleVelocities';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  ND1 = Ncell_p;
  ND3 = Ncell_p^3;
  h5create(foutname,datasetname,[ND3 3]);
  h5write(foutname,datasetname,vc1(:), [1 1], [ND3 1]);
end
clear vc1  %% save memory

%% ------------- vc2 ----------------------
disp('----- Calculating CDM velocity y -----');
vc2(:,:,:)          = -i*af*k2_3D_p./ksq_p.*Thc(:,:,:);
vc2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc2                 = real(ifftn(ifftshift(vc2)));

Vc2 = reshape(vc2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity y more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc2 = padarray(vc2, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc2 = interpn(vc2, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
vc2 = vc2 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vc2'], 'w');
  fwrite(fout, vc2, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,vc2(:), [1 2], [ND3 1]);
end
clear vc2  %% save memory

%% ------------- vc3 ----------------------
disp('----- Calculating CDM velocity z -----');
vc3(:,:,:)          = -i*af*k3_3D_p./ksq_p.*Thc(:,:,:);
vc3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc3                 = real(ifftn(ifftshift(vc3)));

Vc3 = reshape(vc3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity y more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc3 = padarray(vc3, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc3 = interpn(vc3, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
vc3 = vc3 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vc3'], 'w');
  fwrite(fout, vc3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,vc3(:), [1 3], [ND3 1]);
  %% Finalize writing Particle Velocities in hdf5
  topgriddims = -99999*ones(1,3);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(3));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(1));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(ND3));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims-1));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end
clear vc3  %% save memory

%% ------------- velocity divergence ---------------------
Thc = real(ifftn(ifftshift(Thc)));  

ZThc = reshape(Thc(:,:,1),Nmode_p,Nmode_p); %% for figure
clear Thc  %% save memory
%% =========== CDM velocity ==================================== end



%% =========== baryon density ================================== begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
disp('----- Interpolating transfer function -----');
db  = interp2(muext,log(ksampletab), deltasb,  costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
db = rand_real_norm(db,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%%%% If SPH particle is used, one can here get the particle positions 
%%%% just the way CDM positions are calculated.
if baryonparticleflag
  %% ------------- bpos1 ----------------------
  disp('----- Calculating baryon position x -----');
  Psi1                 = i*k1_3D_p./ksq_p.*db;
  Psi1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
  Psi1                 = real(ifftn(ifftshift(Psi1)));
  fout = fopen([ICsubdir '/bpos1'], 'w');
  fwrite(fout, mod((Psi1 + (k1_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  fclose(fout);
  xbar_plane    =   Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
  xbar_ex_plane = 5*Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in baryon position
  if particlevelocity_accuracyflag  
    Psi1 = mod((Psi1 + (k1_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
  else
    clear Psi1  %% save memory
  end
  
  %% ------------- bpos2 ----------------------
  disp('----- Calculating baryon position y -----');
  Psi2                 = i*k2_3D_p./ksq_p.*db;
  Psi2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
  Psi2                 = real(ifftn(ifftshift(Psi2)));
  fout = fopen([ICsubdir '/bpos2'], 'w');
  fwrite(fout, mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  fclose(fout);
  ybar_plane    =   Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
  ybar_ex_plane = 5*Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in CDM position
  if particlevelocity_accuracyflag
    Psi2 = mod((Psi2 + (k2_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
  else
    clear Psi2  %% save memory
  end
  
  %% ------------- bpos3 ----------------------
  disp('----- Calculating baryon position z -----');
  Psi3                 = i*k3_3D_p./ksq_p.*db;
  Psi3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
  Psi3                 = real(ifftn(ifftshift(Psi3)));
  fout = fopen([ICsubdir '/bpos3'], 'w');
  fwrite(fout, mod((Psi3 + (k3_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  fclose(fout);
  if particlevelocity_accuracyflag
    Psi3 = mod((Psi3 + (k3_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
  else
    clear Psi3  %% save memory
  end
end
%% ------------- density ---------------------
disp('----- Calculating baryon density -----');
db = real(ifftn(ifftshift(db)));  

Zb    = reshape(db(:,:,1),Nmode_p,Nmode_p);

%% enzo baryon density output is the following:
%%  db_enzo    = (db+1)*fb

if enzo_bin_flag
  fout = fopen([ICsubdir '/db'], 'w');
  fwrite(fout, (db+1)*fb, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Write Grid Density
  datname     = 'GridDensity';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  h5create(foutname,datasetname,[ND1 ND1 ND1 1]);
  
  griddims    = ND1*ones(1,3);
  topgriddims = ND1*ones(1,3);
  h5write(foutname,datasetname,(db+1)*fb);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(1));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(3));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(griddims));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end

clear db  %% save memory
%% =========== baryon density ================================== end



%% =========== baryon velocity ==================================== begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
disp('----- Interpolating transfer function -----');
Thb = interp2(muext,log(ksampletab), deltasThb, costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
Thb = rand_real_norm(Thb,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vb1 ----------------------
disp('----- Calculating baryon velocity x -----');
vb1(:,:,:)          = -i*af*k1_3D_p./ksq_p.*Thb(:,:,:);
vb1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb1                 = real(ifftn(ifftshift(vb1)));
%% add streaming velocity (V_cb = Vc - Vb)
vb1                 = vb1 - V_cb_1_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--1--**)
sp_Etot_enzo = 1/2*(vb1*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb1 = reshape(vb1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

%% into enzo velocity unit
vb1 = vb1 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb1'], 'w');
  fwrite(fout, vb1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Prepare to write Baryon Grid Velocities in hdf5
  datname     = 'GridVelocities';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  h5create(foutname,datasetname,[ND1 ND1 ND1 3]);
  h5write(foutname,datasetname, vb1, [1 1 1 1], [ND1 ND1 ND1 1]);
end

%% If SPH in mind:
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity x more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb1 = padarray(vb1, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb1 = interpn(vb1, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
  fout = fopen([ICsubdir '/vpb1'], 'w');
  fwrite(fout, vb1, 'double');
  fclose(fout);
end

clear vb1  %% save memory

%% ------------- vb2 ----------------------
disp('----- Calculating baryon velocity y -----');
vb2(:,:,:)          = -i*af*k2_3D_p./ksq_p.*Thb(:,:,:);
vb2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb2                 = real(ifftn(ifftshift(vb2)));
%% add streaming velocity (V_cb = Vc - Vb)
vb2                 = vb2 - V_cb_2_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--2--**)
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb2*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb2 = reshape(vb2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

%% into enzo velocity unit
vb2 = vb2 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb2'], 'w');
  fwrite(fout, vb2, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname, vb2, [1 1 1 2], [ND1 ND1 ND1 1]);
end

if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity y more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb2 = padarray(vb2, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb2 = interpn(vb2, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
  fout = fopen([ICsubdir '/vpb2'], 'w');
  fwrite(fout, vb2, 'double');
  fclose(fout);
end

clear vb2  %% save memory

%% ------------- vb3 ----------------------
disp('----- Calculating baryon velocity z -----');
vb3(:,:,:)          = -i*af*k3_3D_p./ksq_p.*Thb(:,:,:);
vb3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb3                 = real(ifftn(ifftshift(vb3)));
%% add streaming velocity (V_cb = Vc - Vb)
vb3                 = vb3 - V_cb_3_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--3--**)
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb3*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb3 = reshape(vb3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure

%% into enzo velocity unit
vb3 = vb3 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb3'], 'w');
  fwrite(fout, vb3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname, vb3, [1 1 1 3], [ND1 ND1 ND1 1]);
  %% Finalize writing Grid Velocities in hdf5
  griddims    = ND1*ones(1,3);
  topgriddims = ND1*ones(1,3);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(3));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(3));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(griddims));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end

if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity z more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb3 = padarray(vb3, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb3 = interpn(vb3, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
  fout = fopen([ICsubdir '/vpb3'], 'w');
  fwrite(fout, vb3, 'double');
  fclose(fout);
end

clear vb3  %% save memory

%% ------------- velocity divergence ---------------------
Thb = real(ifftn(ifftshift(Thb)));  

ZThb = reshape(Thb(:,:,1),Nmode_p,Nmode_p); %% for figure
clear Thb  %% save memory
%% =========== baryon velocity ==================================== end


%% =========== baryon temperature, energies ======================= begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
disp('----- Interpolating transfer function -----');
dT  = interp2(muext,log(ksampletab), deltasT,  costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
dT = rand_real_norm(dT,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- temperature, energies  ---------------------
disp('----- Calculating baryon temperature and energy -----');
dT = real(ifftn(ifftshift(dT)));  


%% Mean IGM temperature fit from Tseliakhovich & Hirata (2010)
aa1  = 1/119;
aa2  = 1/115;
Tz   = TCMB0/af /(1+af/aa1/(1+(aa2/af)^1.5));  %% in K, global average temperature
Tz   = Tz*(1+DT3D_azend(ic,jc,kc)); %% local(cell) average temperature

mmw = 1.2195; %% mean molecular weight with X=0.76, Y=0.24, neutral.

%% specific thermal energy for monatomic gas (H+He), in units of VelocityUnits^2 : sp_Eth_enzo
%%Tcell = (dT+1)*Tz;  %% in K
%%sp_Eth_enzo  = 3/2*kb*Tcell /(mmw*mH) /VelocityUnits^2;  %% see Enzo paper(2014) eq. 7.
if enzo_bin_flag
  fout = fopen([ICsubdir '/etherm'], 'w');
  fwrite(fout, 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Write Baryon Thermal Energy
  datname     = 'GasThermalSpecEnergy';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  h5create(foutname,datasetname,[ND1 ND1 ND1 1]);

  griddims    = ND1*ones(1,3);
  topgriddims = ND1*ones(1,3);
  h5write(foutname,datasetname,3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(1));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(3));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(griddims));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end

%% in erg (thermal energy per baryon)
Zeth = reshape(3/2*kb*(dT(:,:,1)+1)*Tz/(mmw*mH),Nmode_p,Nmode_p); %% for figure

%% in K
Ztemp = reshape((dT(:,:,1)+1)*Tz,Nmode_p,Nmode_p); %% for figure


%% specific total energy for monatomic gas (H+He), in units of VelocityUnits^2 :   sp_Etot_enzo
%% Currently sp_Etot_enzo does not include magnetic contribution, but in principle it should.
%% memory-saving way of calculating sp_Etot_enzo (**--4--**)
sp_Etot_enzo = sp_Etot_enzo + 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2;
  
if enzo_bin_flag
  fout = fopen([ICsubdir '/etot'], 'w');
  fwrite(fout, sp_Etot_enzo, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  %% Write Baryon Total Energy
  datname     = 'GasTotalSpecEnergy';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); %% In case file already exists, delete the file
  h5create(foutname,datasetname,[ND1 ND1 ND1 1]);

  griddims    = ND1*ones(1,3);
  topgriddims = ND1*ones(1,3);
  h5write(foutname,datasetname,sp_Etot_enzo);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(1));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(3));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(griddims));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end

%% in erg (total energy per baryon)
Zetot = reshape(sp_Etot_enzo(:,:,1)*VelocityUnits^2,Nmode_p,Nmode_p); %% for figure

%% If SPH in mind: (this if-end statement should be placed after other dT related
%% grid quantities are all calculated. So developers, do not change this location.)
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle thermal energy more accurately than 1LPT **');
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  dT = padarray(dT, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  dT = interpn(dT, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
  fout = fopen([ICsubdir '/eptherm'], 'w');
  fwrite(fout, 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2, 'double');
  fclose(fout);
end

clear dT sp_Etot_enzo Psi1 Psi2 Psi3
