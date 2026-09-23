%% This script generates initial condition data for enzo and dumps them.
%% For ICs for other simulation codes, this is the script to start from.
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
disp('----- Interpolating transfer function -----');
dc  = zeros(Nmode_p,Nmode_p,Nmode_p);
if ~memory_save
  dc  = interp2(muext,log(ksampletab), deltasc,  costh_k_V, 0.5*log(ksq_p),interp2opt);  %% dc still k-space values here.
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    %% k3_3D_p_chunk is defined at the bottom chunk, so need to add a shift in iteration
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      dc(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasc,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2),interp2opt);  %% dc still k-space values
    else %% the last chunk
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      dc(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasc,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2),interp2opt);  %% dc still k-space values
    end    
  end
end

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
dc = rand_real_norm(dc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- cpos1 ----------------------
disp('----- Calculating CDM position x -----');
if ~memory_save
  Psi1                 = 1i*k1_3D_p./ksq_p.*dc;
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      Psi1(:,:,kkstart:kkend) = 1i*k1_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*dc(:,:,kkstart:kkend);
    else %% the last chunk
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      Psi1(:,:,kkstart:kkend) = 1i*k1_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*dc(:,:,kkstart:kkend);
    end
  end
end
Psi1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi1                 = real(ifftn(ifftshift(Psi1)));

if ~memory_save
  xCDM_plane    =   Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
  xCDM_ex_plane = 5*Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
else
  xCDM_plane    =   Psi1(:,:,1) + k1_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
  xCDM_ex_plane = 5*Psi1(:,:,1) + k1_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
end

if ~memory_save
  Psi1 = mod((Psi1 + (k1_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    if kkchunk ~= Nchunk
      Psi1(:,:,kkstart:kkend) = mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk            /kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Psi1(:,:,kkstart:kkend) = mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk(:,:,1:Nsub)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    end
  end
end

if enzo_bin_flag
  fout = fopen([ICsubdir '/cpos1'], 'w');
  fwrite(fout, Psi1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  datname     = 'ParticlePositions';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname);
  ND1 = Ncell_p;
  ND3 = Ncell_p^3;
  h5create(foutname,datasetname,[ND3 3]);
  h5write(foutname,datasetname,Psi1(:), [1 1], [ND3 1]);
end

if particlevelocity_accuracyflag
    Psi1 = mod(Psi1*Nmode_p - 0.5, Nmode_p)+1;
else
  clear Psi1  
end

%% ------------- cpos2 ----------------------
disp('----- Calculating CDM position y -----');
if ~memory_save
  Psi2                 = 1i*k2_3D_p./ksq_p.*dc;
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      Psi2(:,:,kkstart:kkend) = 1i*k2_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*dc(:,:,kkstart:kkend);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Psi2(:,:,kkstart:kkend) = 1i*k2_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*dc(:,:,kkstart:kkend);
    end
  end
end
Psi2(Nc_p,Nc_p,Nc_p) = complex(0); 
Psi2                 = real(ifftn(ifftshift(Psi2)));

if ~memory_save
  yCDM_plane    =   Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  yCDM_ex_plane = 5*Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
else
  yCDM_plane    =   Psi2(:,:,1) + k2_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  yCDM_ex_plane = 5*Psi2(:,:,1) + k2_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
end

if ~memory_save
  Psi2 = mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    if kkchunk ~= Nchunk
      Psi2(:,:,kkstart:kkend) = mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk            /kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Psi2(:,:,kkstart:kkend) = mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk(:,:,1:Nsub)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    end
  end
end

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
  clear Psi2  
end

%% ------------- cpos3 ----------------------
disp('----- Calculating CDM position z -----');
if ~memory_save
  Psi3                 = 1i*k3_3D_p./ksq_p.*dc;
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      Psi3(:,:,kkstart:kkend) = 1i*(k3_3D_p_chunk            +kshift)./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*dc(:,:,kkstart:kkend);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Psi3(:,:,kkstart:kkend) = 1i*(k3_3D_p_chunk(:,:,1:Nsub)+kshift)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*dc(:,:,kkstart:kkend);
    end
  end
end
Psi3(Nc_p,Nc_p,Nc_p) = complex(0); 
Psi3                 = real(ifftn(ifftshift(Psi3)));

if ~memory_save
  zCDM_plane    =   Psi3(:,:,1) + k3_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  zCDM_ex_plane = 5*Psi3(:,:,1) + k3_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
else
  zCDM_plane    =   Psi3(:,:,1) + k3_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  zCDM_ex_plane = 5*Psi3(:,:,1) + k3_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
end

if ~memory_save
  Psi3 = mod((Psi3 + (k3_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;  %% k3_3D_p_chunk holds the first chunk only, so shift k3 by kshift for each chunk.
    if kkchunk ~= Nchunk
      Psi3(:,:,kkstart:kkend) = mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk            +kshift)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Psi3(:,:,kkstart:kkend) = mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk(:,:,1:Nsub)+kshift)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1);
    end
  end
end

if enzo_bin_flag
  fout = fopen([ICsubdir '/cpos3'], 'w');
  fwrite(fout, Psi3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,Psi3(:), [1 3], [ND3 1]);
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
  clear Psi3 
end

%% ------------- density ---------------------
dc = real(ifftn(ifftshift(dc))); 
Zc = reshape(dc(:,:,1),Nmode_p,Nmode_p); 
clear dc  

%% =========== CDM velocity ==================================== begin
disp('----- Interpolating transfer function -----');
Thc  = zeros(Nmode_p,Nmode_p,Nmode_p);
if ~memory_save
  Thc  = interp2(muext,log(ksampletab), deltasThc,  costh_k_V, 0.5*log(ksq_p),interp2opt); 
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      Thc(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasThc,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2),interp2opt); 
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      Thc(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasThc,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2),interp2opt); 
    end    
  end
end

%% randomize, apply reality, and normalize
disp('----- Convolving transfer function with random number -----');
Thc = rand_real_norm(Thc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vc1 ----------------------
disp('----- Calculating CDM velocity x -----');
if ~memory_save
  vc1(:,:,:)          = -1i*af*k1_3D_p./ksq_p.*Thc(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vc1(:,:,kkstart:kkend) = -1i*af*k1_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thc(:,:,kkstart:kkend);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      vc1(:,:,kkstart:kkend) = -1i*af*k1_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thc(:,:,kkstart:kkend);
    end
  end
end

vc1(Nc_p,Nc_p,Nc_p) = complex(0); 
vc1                 = real(ifftn(ifftshift(vc1)));
Vc1 = reshape(vc1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity x more accurately than 1LPT **');
  vc1 = padarray(vc1, [1 1 1], 'circular', 'post'); 
  vc_1 = zeros(Nmode_p,Nmode_p,Nmode_p); 
  
  if ~memory_save
    vc_1 = interpn(vc1, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      vc_1(:,:,kkstart:kkend) = interpn(vc1, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vc1 = vc_1;
  clear vc_1;
end

vc1 = vc1 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vc1'], 'w');
  fwrite(fout, vc1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  datname     = 'ParticleVelocities';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); 
  ND1 = Ncell_p;
  ND3 = Ncell_p^3;
  h5create(foutname,datasetname,[ND3 3]);
  h5write(foutname,datasetname,vc1(:), [1 1], [ND3 1]);
end
clear vc1  

%% ------------- vc2 ----------------------
disp('----- Calculating CDM velocity y -----');
if ~memory_save
  vc2(:,:,:)          = -1i*af*k2_3D_p./ksq_p.*Thc(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vc2(:,:,kkstart:kkend) = -1i*af*k2_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thc(:,:,kkstart:kkend);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      vc2(:,:,kkstart:kkend) = -1i*af*k2_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thc(:,:,kkstart:kkend);
    end
  end
end

vc2(Nc_p,Nc_p,Nc_p) = complex(0); 
vc2                 = real(ifftn(ifftshift(vc2)));
Vc2 = reshape(vc2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity y more accurately than 1LPT **');
  vc2 = padarray(vc2, [1 1 1], 'circular', 'post'); 
  vc_2 = zeros(Nmode_p,Nmode_p,Nmode_p); 
  
  if ~memory_save
    vc_2 = interpn(vc2, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      vc_2(:,:,kkstart:kkend) = interpn(vc2, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vc2 = vc_2;
  clear vc_2;
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
clear vc2  

%% ------------- vc3 ----------------------
disp('----- Calculating CDM velocity z -----');
if ~memory_save
  vc3(:,:,:)          = -1i*af*k3_3D_p./ksq_p.*Thc(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vc3(:,:,kkstart:kkend) = -1i*af*(k3_3D_p_chunk            +kshift)./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thc(:,:,kkstart:kkend);
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      vc3(:,:,kkstart:kkend) = -1i*af*(k3_3D_p_chunk(:,:,1:Nsub)+kshift)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thc(:,:,kkstart:kkend);
    end
  end
end

vc3(Nc_p,Nc_p,Nc_p) = complex(0); 
vc3                 = real(ifftn(ifftshift(vc3)));
Vc3 = reshape(vc3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

if particlevelocity_accuracyflag
  disp('******* Calculating CDM velocity z more accurately than 1LPT **');
  vc3 = padarray(vc3, [1 1 1], 'circular', 'post'); 
  vc_3 = zeros(Nmode_p,Nmode_p,Nmode_p); 
  
  if ~memory_save
    vc_3 = interpn(vc3, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      vc_3(:,:,kkstart:kkend) = interpn(vc3, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vc3 = vc_3;
  clear vc_3;
end

vc3 = vc3 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vc3'], 'w');
  fwrite(fout, vc3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname,vc3(:), [1 3], [ND3 1]);
  topgriddims = -99999*ones(1,3);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(3));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(1));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(ND3));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims-1));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end
clear vc3  

%% ------------- velocity divergence ---------------------
Thc = real(ifftn(ifftshift(Thc)));  
ZThc = reshape(Thc(:,:,1),Nmode_p,Nmode_p); 
clear Thc  

%% =========== baryon density ================================== begin
disp('----- Interpolating transfer function -----');
db  = zeros(Nmode_p,Nmode_p,Nmode_p);
if ~memory_save
  db  = interp2(muext,log(ksampletab), deltasb,  costh_k_V, 0.5*log(ksq_p),interp2opt); 
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      db(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasb,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2),interp2opt); 
    else
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      db(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasb,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2),interp2opt); 
    end    
  end
end

disp('----- Convolving transfer function with random number -----');
db = rand_real_norm(db,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

if baryonparticleflag
  disp('----- Calculating baryon position x -----');
  if ~memory_save
    Psi1                 = 1i*k1_3D_p./ksq_p.*db;
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
      if kkchunk ~= Nchunk
        Psi1(:,:,kkstart:kkend) = 1i*k1_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*db(:,:,kkstart:kkend);
      else
         if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        Psi1(:,:,kkstart:kkend) = 1i*k1_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*db(:,:,kkstart:kkend);
      end
    end
  end
  Psi1(Nc_p,Nc_p,Nc_p) = complex(0); 
  Psi1                 = real(ifftn(ifftshift(Psi1)));
  fout = fopen([ICsubdir '/bpos1'], 'w');
  if ~memory_save
    fwrite(fout, mod((Psi1 + (k1_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      if kkchunk ~= Nchunk
        fwrite(fout, mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk            /kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      else
        if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        fwrite(fout, mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk(:,:,1:Nsub)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      end
    end
  end
  fclose(fout);
  if ~memory_save
    xbar_plane    =   Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
    xbar_ex_plane = 5*Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  else
    xbar_plane    =   Psi1(:,:,1) + k1_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
    xbar_ex_plane = 5*Psi1(:,:,1) + k1_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  end

  if particlevelocity_accuracyflag
      if ~memory_save
      Psi1 = mod((Psi1 + (k1_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
    else
      for kkchunk=1:Nchunk
        kkstart = 1 + (kkchunk-1)*Nsub_chunk;
        kkend   = kkchunk*Nsub_chunk;
        if kkchunk == Nchunk; kkend=Nmode_p; end
        if kkchunk ~= Nchunk
          Psi1(:,:,kkstart:kkend) = mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk            /kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        else 
          if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
          Psi1(:,:,kkstart:kkend) = mod((Psi1(:,:,kkstart:kkend) + (k1_3D_p_chunk(:,:,1:Nsub)/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        end
      end
    end
  else
    clear Psi1  
  end

  disp('----- Calculating baryon position y -----');
  if ~memory_save
    Psi2                 = 1i*k2_3D_p./ksq_p.*db;
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
      if kkchunk ~= Nchunk
        Psi2(:,:,kkstart:kkend) = 1i*k2_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*db(:,:,kkstart:kkend);
      else
        if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        Psi2(:,:,kkstart:kkend) = 1i*k2_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*db(:,:,kkstart:kkend);
      end
    end
  end
  Psi2(Nc_p,Nc_p,Nc_p) = complex(0); 
  Psi2                 = real(ifftn(ifftshift(Psi2)));
  fout = fopen([ICsubdir '/bpos2'], 'w');
  if ~memory_save
    fwrite(fout, mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      if kkchunk ~= Nchunk
        fwrite(fout, mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk            /kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      else
        if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        fwrite(fout, mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk(:,:,1:Nsub)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      end
    end
  end
  fclose(fout);
  if ~memory_save
    ybar_plane    =   Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
    ybar_ex_plane = 5*Psi2(:,:,1) + k2_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  else
    ybar_plane    =   Psi2(:,:,1) + k2_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
    ybar_ex_plane = 5*Psi2(:,:,1) + k2_3D_p_chunk(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; 
  end

  if particlevelocity_accuracyflag
    if ~memory_save
      Psi2 = mod((Psi2 + (k2_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
    else
      for kkchunk=1:Nchunk
        kkstart = 1 + (kkchunk-1)*Nsub_chunk;
        kkend   = kkchunk*Nsub_chunk;
        if kkchunk == Nchunk; kkend=Nmode_p; end
        if kkchunk ~= Nchunk
          Psi2(:,:,kkstart:kkend) = mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk             /kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        else 
          if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
          Psi2(:,:,kkstart:kkend) = mod((Psi2(:,:,kkstart:kkend) + (k2_3D_p_chunk(:,:,1:Nsub)/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        end
      end
    end
  else
    clear Psi2 
  end

  disp('----- Calculating baryon position z -----');
  if ~memory_save
    Psi3                 = 1i*k3_3D_p./ksq_p.*db;
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
      if kkchunk ~= Nchunk
        Psi3(:,:,kkstart:kkend) = 1i*(k3_3D_p_chunk            +kshift)./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*db(:,:,kkstart:kkend);
      else 
        if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        Psi3(:,:,kkstart:kkend) = 1i*(k3_3D_p_chunk(:,:,1:Nsub)+kshift)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*db(:,:,kkstart:kkend);
      end
    end
  end
  Psi3(Nc_p,Nc_p,Nc_p) = complex(0); 
  Psi3                 = real(ifftn(ifftshift(Psi3)));
  fout = fopen([ICsubdir '/bpos3'], 'w');
  if ~memory_save
    fwrite(fout, mod((Psi3 + (k3_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  else
    for kkchunk=1:Nchunk
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk; kkend=Nmode_p; end
      kshift = kunit_p*(kkchunk-1)*Nsub_chunk;  %% CRITICAL FIX: kshift recalculation
      if kkchunk ~= Nchunk
        fwrite(fout, mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk            +kshift)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      else 
        if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
        fwrite(fout, mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk(:,:,1:Nsub)+kshift)/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
      end
    end
  end
  fclose(fout);

  if particlevelocity_accuracyflag
    if ~memory_save
      Psi3 = mod((Psi3 + (k3_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
    else
      for kkchunk=1:Nchunk
        kkstart = 1 + (kkchunk-1)*Nsub_chunk;
        kkend   = kkchunk*Nsub_chunk;
        if kkchunk == Nchunk
          kkend=Nmode_p;
        end
        kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
        if kkchunk ~= Nchunk
          Psi3(:,:,kkstart:kkend) = mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk            +kshift)/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        else 
          if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
          Psi3(:,:,kkstart:kkend) = mod((Psi3(:,:,kkstart:kkend) + ((k3_3D_p_chunk(:,:,1:Nsub)+kshift)/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
        end
      end
    end
  else
    clear Psi3  
  end
end

disp('----- Calculating baryon density -----');
db = real(ifftn(ifftshift(db)));  
Zb    = reshape(db(:,:,1),Nmode_p,Nmode_p); 

%% --- SAFETY GUARD: Use local baryon fraction for internal grid consistency ---
fb_l = (1+cellspec_azend(idxcc,5))*fb / ((1+cellspec_azend(idxcc,4))*fc + (1+cellspec_azend(idxcc,5))*fb);

if enzo_bin_flag
  fout = fopen([ICsubdir '/db'], 'w');
  fwrite(fout, (db+1)*fb_l, 'double');
  fclose(fout);
end

if enzo_HDF5_flag
  datname     = 'GridDensity';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); 

  h5create(foutname,datasetname,[ND1 ND1 ND1 1]);
  griddims    = ND1*ones(1,3);
  topgriddims = ND1*ones(1,3);
  h5write(foutname,datasetname,(db+1)*fb_l);
  h5writeatt(foutname,datasetname,'Component_Rank',int64(1));
  h5writeatt(foutname,datasetname,'Component_Size',int64(ND3));
  h5writeatt(foutname,datasetname,'Rank',          int64(3));
  h5writeatt(foutname,datasetname,'Dimensions',    int64(griddims));
  h5writeatt(foutname,datasetname,'TopGridDims',   int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridEnd',    int64(topgriddims));
  h5writeatt(foutname,datasetname,'TopGridStart',  int64(zeros(1,3)));
end
clear db  


%% =========== baryon velocity ==================================== begin
disp('----- Interpolating transfer function -----');
Thb  = zeros(Nmode_p,Nmode_p,Nmode_p);
if ~memory_save
  Thb  = interp2(muext,log(ksampletab), deltasThb,  costh_k_V, 0.5*log(ksq_p),interp2opt); 
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      Thb(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasThb,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2),interp2opt); 
    else 
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      Thb(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasThb,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2),interp2opt); 
    end    
  end
end

disp('----- Convolving transfer function with random number -----');
Thb = rand_real_norm(Thb,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vb1 ----------------------
disp('----- Calculating baryon velocity x -----');
if ~memory_save
  vb1(:,:,:)          = -1i*af*k1_3D_p./ksq_p.*Thb(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vb1(:,:,kkstart:kkend) = -1i*af*k1_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thb(:,:,kkstart:kkend);
    else 
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      vb1(:,:,kkstart:kkend) = -1i*af*k1_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thb(:,:,kkstart:kkend);
    end
  end
end

vb1(Nc_p,Nc_p,Nc_p) = complex(0); 
vb1                 = real(ifftn(ifftshift(vb1)));

vb1                 = vb1 - V_cb_1_azend(ic,jc,kc); 
sp_Etot_enzo = 1/2*(vb1*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb1 = reshape(vb1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

vb1 = vb1 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb1'], 'w');
  fwrite(fout, vb1, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  datname     = 'GridVelocities';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); 
  h5create(foutname,datasetname,[ND1 ND1 ND1 3]);
  h5write(foutname,datasetname, vb1, [1 1 1 1], [ND1 ND1 ND1 1]);
end

%% If SPH in mind:
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity x more accurately than 1LPT **');
  vb1 = padarray(vb1, [1 1 1], 'circular', 'post'); 
  vb_1 = zeros(Nmode_p,Nmode_p,Nmode_p); %% temporary variable. Allocated only when used.
  if ~memory_save
    vb_1 = interpn(vb1, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk kkend=Nmode_p; end
      vb_1(:,:,kkstart:kkend) = interpn(vb1, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vb1 = vb_1;
  clear vb_1;
  fout = fopen([ICsubdir '/vpb1'], 'w');
  fwrite(fout, vb1, 'double');
  fclose(fout);
end
clear vb1  

%% ------------- vb2 ----------------------
disp('----- Calculating baryon velocity y -----');
if ~memory_save
  vb2(:,:,:)          = -1i*af*k2_3D_p./ksq_p.*Thb(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vb2(:,:,kkstart:kkend) = -1i*af*k2_3D_p_chunk            ./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thb(:,:,kkstart:kkend);
    else 
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      vb2(:,:,kkstart:kkend) = -1i*af*k2_3D_p_chunk(:,:,1:Nsub)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thb(:,:,kkstart:kkend);
    end
  end
end

vb2(Nc_p,Nc_p,Nc_p) = complex(0); 
vb2                 = real(ifftn(ifftshift(vb2)));

vb2                 = vb2 - V_cb_2_azend(ic,jc,kc); 
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb2*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb2 = reshape(vb2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

vb2 = vb2 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb2'], 'w');
  fwrite(fout, vb2, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname, vb2, [1 1 1 2], [ND1 ND1 ND1 1]);
end

%% If SPH in mind:
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity y more accurately than 1LPT **');
  vb2 = padarray(vb2, [1 1 1], 'circular', 'post'); 
  vb_2 = zeros(Nmode_p,Nmode_p,Nmode_p); %% CRITICAL FIX: Memory allocation inside block
  if ~memory_save
    vb_2 = interpn(vb2, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk kkend=Nmode_p; end
      vb_2(:,:,kkstart:kkend) = interpn(vb2, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vb2 = vb_2;
  clear vb_2;
  fout = fopen([ICsubdir '/vpb2'], 'w');
  fwrite(fout, vb2, 'double');
  fclose(fout);
end
clear vb2  

%% ------------- vb3 ----------------------
disp('----- Calculating baryon velocity z -----');
if ~memory_save
  vb3(:,:,:)          = -1i*af*k3_3D_p./ksq_p.*Thb(:,:,:);
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend=Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      vb3(:,:,kkstart:kkend) = -1i*af*(k3_3D_p_chunk            +kshift)./(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2).*Thb(:,:,kkstart:kkend);
    else 
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      vb3(:,:,kkstart:kkend) = -1i*af*(k3_3D_p_chunk(:,:,1:Nsub)+kshift)./(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2).*Thb(:,:,kkstart:kkend);
    end
  end
end

vb3(Nc_p,Nc_p,Nc_p) = complex(0); 
vb3                 = real(ifftn(ifftshift(vb3)));

vb3                 = vb3 - V_cb_3_azend(ic,jc,kc); 
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb3*MpcMyr_2_kms*1e5/VelocityUnits).^2; 

Vb3 = reshape(vb3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); 

vb3 = vb3 * MpcMyr_2_kms * 1e5 /VelocityUnits;

if enzo_bin_flag
  fout = fopen([ICsubdir '/vb3'], 'w');
  fwrite(fout, vb3, 'double');
  fclose(fout);
end
if enzo_HDF5_flag
  h5write(foutname,datasetname, vb3, [1 1 1 3], [ND1 ND1 ND1 1]);
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

%% If SPH in mind:
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle velocity z more accurately than 1LPT **');
  vb3 = padarray(vb3, [1 1 1], 'circular', 'post'); 
  vb_3 = zeros(Nmode_p,Nmode_p,Nmode_p); %% CRITICAL FIX: Memory allocation inside block
  if ~memory_save
    vb_3 = interpn(vb3, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk kkend=Nmode_p; end
      vb_3(:,:,kkstart:kkend) = interpn(vb3, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  vb3 = vb_3;
  clear vb_3;
  fout = fopen([ICsubdir '/vpb3'], 'w');
  fwrite(fout, vb3, 'double');
  fclose(fout);
end
clear vb3  

Thb = real(ifftn(ifftshift(Thb)));  
ZThb = reshape(Thb(:,:,1),Nmode_p,Nmode_p); 
clear Thb  

%% =========== baryon temperature, energies ======================= begin
disp('----- Interpolating transfer function -----');
dT  = zeros(Nmode_p,Nmode_p,Nmode_p);
if ~memory_save
  dT  = interp2(muext,log(ksampletab), deltasT,  costh_k_V, 0.5*log(ksq_p),interp2opt); 
else
  for kkchunk=1:Nchunk
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk; kkend=Nmode_p; end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      dT(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasT,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2),interp2opt); 
    else 
      if Nsub_r == 0; Nsub = Nsub_chunk; else; Nsub = Nsub_r; end
      dT(:,:,kkstart:kkend) = interp2(muext,log(ksampletab), deltasT,  costh_k_V(:,:,kkstart:kkend), 0.5*log(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2),interp2opt); 
    end    
  end
end

disp('----- Convolving transfer function with random number -----');
dT = rand_real_norm(dT,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

disp('----- Calculating baryon temperature and energy -----');
dT = real(ifftn(ifftshift(dT)));  

aa1  = 1/119;
aa2  = 1/115;
Tz   = TCMB0/af /(1+af/aa1/(1+(aa2/af)^1.5));  
Tz   = Tz*(1+DT3D_azend(ic,jc,kc)); 
mmw = 1.2195; 

if enzo_bin_flag
  fout = fopen([ICsubdir '/etherm'], 'w');
  fwrite(fout, 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2, 'double');
  fclose(fout);
end

if enzo_HDF5_flag
  datname     = 'GasThermalSpecEnergy';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); 
  
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

Zeth = reshape(3/2*kb*(dT(:,:,1)+1)*Tz/(mmw*mH),Nmode_p,Nmode_p); 
Ztemp = reshape((dT(:,:,1)+1)*Tz,Nmode_p,Nmode_p); 

sp_Etot_enzo = sp_Etot_enzo + 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2;

if enzo_bin_flag
  fout = fopen([ICsubdir '/etot'], 'w');
  fwrite(fout, sp_Etot_enzo, 'double');
  fclose(fout);
end

if enzo_HDF5_flag
  datname     = 'GasTotalSpecEnergy';
  foutname    = [ICsubdir '/' datname];
  datasetname = ['/' datname];
  delete(foutname); 
  
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

Zetot = reshape(sp_Etot_enzo(:,:,1)*VelocityUnits^2,Nmode_p,Nmode_p); 

%% If SPH in mind: (this if-end statement should be placed after other dT related
%% grid quantities are all calculated. So developers, do not change this location.)
if (particlevelocity_accuracyflag & baryonparticleflag)
  disp('******* Calculating baryon particle thermal energy more accurately than 1LPT **');
  dT = padarray(dT, [1 1 1], 'circular', 'post'); 
  dT_ = zeros(Nmode_p,Nmode_p,Nmode_p); %% CRITICAL FIX: Memory allocation inside block
  if ~memory_save
    dT_ = interpn(dT, Psi1, Psi2, Psi3, interpnopt);
  else
    for kkchunk=1:Nchunk  
      kkstart = 1 + (kkchunk-1)*Nsub_chunk;
      kkend   = kkchunk*Nsub_chunk;
      if kkchunk == Nchunk kkend=Nmode_p; end
      dT_(:,:,kkstart:kkend) = interpn(dT, Psi1(:,:,kkstart:kkend), Psi2(:,:,kkstart:kkend), Psi3(:,:,kkstart:kkend), interpnopt);
    end
  end
  dT = dT_;
  clear dT_;
  fout = fopen([ICsubdir '/eptherm'], 'w');
  fwrite(fout, 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2, 'double');
  fclose(fout);
end

clear dT sp_Etot_enzo Psi1 Psi2 Psi3
