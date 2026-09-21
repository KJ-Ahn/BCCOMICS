%% BCCOMICS_CR: Reads in fluctuation (curvature*TF) made by bccomics_CR_setup.m.
%%              Generates 3D fields of small-scale (inside-a-patch)
%%              perturbations into bare binary files.
%%              Modified to bypass patch selection and load Constrained Realization (CR) target.

clear;  %% Clears the memory and have a fresh start!
more off; %% enables to see progress
returnflag = false; %% main program need to stop when script stops.
%% Start recording log
diary on;

%% Detect which is running: octave or matlab?
if (exist('OCTAVE_VERSION','builtin'))
  matlabflag=false;
  disp('----------------run on OCTAVE----------------');
else
  matlabflag=true;
  disp('----------------run on MATLAB----------------');
end

%% Read in constants in cgs unit and conversion factors.
Consts_Conversions;  %%==== script ==================
%% Read in parameters
run('params.m');  %%==== script ==================

%% Check function availability and provide cure
Check_functions;  %%==== script ==================
if returnflag
  clear;
  return;
end

%% Read in cosmology
run(Cosmology);  %%==== script ==================
fb = ombh2/(ombh2+omch2); %% baryon/matter fraction
fc = omch2/(ombh2+omch2); %% CDM/matter fraction

%% IC-for-which-simulation-code flags: default to false, and
%% overwrite them with params_patch.m.
enzo_bin_flag   = false;
enzo_HDF5_flag  = false;
gadget_bin_flag = false;

%% Read in parameters for initial condition
run('params_patch.m');  %%==== script ==================

%% Requires mod(Ncell_p,4)=0 to properly use existing random seed.
if (mod(Ncell_p,4)~=0)
  disp('Choose a number which is multiple of 4 for Ncell_p');
  clear;
  return;
end

%% Setting resolution etc.
%% Read in parameters for initial condition
patch_init;  %%==== script ==================

interp2opt = 'cubic';
interpnopt = 'linear';

%% k1 component on each (k1,k2,k3) point, etc.
if Nmode_p<=32
  memory_save = false;
end

%% number of chunks to iterate for memory_save
if memory_save
  Nsub_chunk = 8;  %% # of slices to treat in one chunk
  Nchunk = ceil(Nmode_p/Nsub_chunk); 
  %% check for residual slices after chunkening
  Nsub_r = mod(Nmode_p,Nsub_chunk); %% # of residual slices
end

if ~memory_save
  [k1_3D_p, k2_3D_p, k3_3D_p] = ndgrid(-Nhalf_p:Nhalf_p-1);
  k1_3D_p = kunit_p * k1_3D_p;
  k2_3D_p = kunit_p * k2_3D_p;
  k3_3D_p = kunit_p * k3_3D_p;
  ksq_p   = k1_3D_p.^2 +k2_3D_p.^2 +k3_3D_p.^2; %% |k|^2
else
  [k1_2D_p, k2_2D_p] = ndgrid(-Nhalf_p:Nhalf_p-1);
  k1_2D_p = kunit_p * k1_2D_p;
  k2_2D_p = kunit_p * k2_2D_p;
  k1_3D_p_chunk = k1_2D_p;
  k2_3D_p_chunk = k2_2D_p;
  for kkk=2:Nsub_chunk
    k1_3D_p_chunk = cat(3, k1_3D_p_chunk, k1_2D_p);
    k2_3D_p_chunk = cat(3, k2_3D_p_chunk, k2_2D_p);
  end
  onesk1k2_2D = ones(Nmode_p,Nmode_p);  
  k3_3D_p_chunk = kunit_p*(-Nhalf_p)*onesk1k2_2D;  
  for kkk=-Nhalf_p+1:-Nhalf_p+(Nsub_chunk-1)
    k3_3D_p_chunk = cat(3, k3_3D_p_chunk, kunit_p*kkk*onesk1k2_2D);
  end
  k1k2sq_p_chunk = k1_3D_p_chunk.^2 +k2_3D_p_chunk.^2; 
end

%% read in mu info
mu  = load([setupdir '/mu.dat']);
dmu = mu(2)-mu(1);
Nmu = length(mu);

%% choose patch to generate initial condition on (CR target bypass)
Load_CR_patch;  %%==== script ==================
if returnflag
  clear;
  return;
end

%% open transfer function file for given patch
ic   = cellspec(idxcc,1);
jc   = cellspec(idxcc,2);
kc   = cellspec(idxcc,3);
strD = [setupdir '/deltas/Deltas_1Dmu_ic' num2str(ic) '_jc' num2str(jc) '_kc' num2str(kc) '-muhalf.matbin'];

if matlabflag
  load(strD, '-mat', 'ksampletab', 'deltasc', 'deltasb', 'deltasThc', 'deltasThb', 'deltasT');
else
  load('-mat-binary', strD, 'ksampletab', 'deltasc', 'deltasb', 'deltasThc', 'deltasThb', 'deltasT');
end

%% Generate initial condition directory
if ~exist(ICdir)
  mkdir(ICdir);
end

Lbox_p_inMpch = Lbox_p*h;  %% enzo uses 'ComovingBoxSize' in units of Mpc/h

%% read in z=zi=1000 statistics for stdDm
fin=fopen([setupdir '/stats_zi.dat']);
fgets(fin); %% skip a line
fgets(fin); %% skip another line
statszi = fscanf(fin, '%e %e %e %e %e %e %e %e %e %e');
fclose(fin);
stdDm_zi = statszi(10);

%% get Dm(in sigma), Vcb, and relative direction for directory name
vcb_val   = cellspec(idxcc, 11);
dm_val    = fc * cellspec(idxcc, 4) + fb * cellspec(idxcc, 5);
sigma_val = dm_val / stdDm_zi;

%% calculate relative direction ratios (e.g., [1, 1, 1])
vdir_raw = -cellspec(idxcc, 8:10);
if vcb_val > 1e-6
    nonzero_v = abs(vdir_raw(abs(vdir_raw) > 1e-6));
    vdir_rel  = round((vdir_raw / min(nonzero_v)) * 10) / 10;
else
    vdir_rel  = [0, 0, 0];
end

dir_str = ['_dir' num2str(vdir_rel(1), '%g') '_' num2str(vdir_rel(2), '%g') '_' num2str(vdir_rel(3), '%g')];

ICsubdir = [ICdir '/' num2str(Lbox_p_inMpch,'%.2f') 'Mpch_' num2str(Ncell_p) '_Dm' num2str(sigma_val, '%.1f') 's_Vcb' num2str(vcb_val, '%.1f') dir_str];
if ~exist(ICsubdir)
  mkdir(ICsubdir); 
end

%% Copy essential files to ICsubdir. This helps consistency in record tracking.
if exist('params.m', 'file'), copyfile('params.m', ICsubdir); end
if exist('params_patch.m', 'file'), copyfile('params_patch.m', ICsubdir); end
if exist(Cosmology, 'file'), copyfile(Cosmology, ICsubdir); end
if exist([setupdir '/zz.dat'], 'file'), copyfile([setupdir '/zz.dat'], ICsubdir); end
if exist([setupdir '/stats_zend.dat'], 'file'), copyfile([setupdir '/stats_zend.dat'], ICsubdir); end

iccdat = cellspec_azend(idxcc,:);
fout   = fopen([ICsubdir '/icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat'],'w');
fprintf(fout, '%4i %4i %4i %e %e %e %e %e %e %e %e %e\n', iccdat');
fclose(fout);

zf = zzend;  
af = 1/(1+zf);  

%% prepare for initial conditions for enzo (set units)
Prepare_enzoIC;  %%==== script ==================

%% Set gaussian random seed
Set_gaussrand;  %%==== script ==================

%% record random seed if wanted
if recordseedflag
  fileNseed = [ICsubdir '/subgaussseed' num2str(Nmode_p) '.matbin'];
  disp(['--- Seed is being recorded under ' ICsubdir ' ---']);
  if matlabflag
    save(fileNseed, 'randamp', 'randphs', '-v6');
  else
    save('-mat-binary', fileNseed, 'randamp', 'randphs');
  end

  if matlabflag
    msg   = dir(fileNseed);
    fsize = msg.bytes; 
  else
    msg   = lstat(fileNseed);
    fsize = msg.size; 
  end  
  if (fsize < Nmode_p*Nmode_p*Nc_p*2*8)  
    delete(fileNseed);
    fileNseed_1 = [ICsubdir '/subgaussseed' num2str(Nmode_p) '.bin'];
    ffout = fopen(fileNseed_1,'w');
    fwrite(ffout, randamp, 'double');
    fwrite(ffout, randphs, 'double');
    fclose(ffout);
  end
end

%% Switching Real and Imag is done easily by 1i*conj(complex_number).
deltasc  (:,Nmu:2*Nmu-1) = deltasc  (:,:);
deltasb  (:,Nmu:2*Nmu-1) = deltasb  (:,:);
deltasThc(:,Nmu:2*Nmu-1) = deltasThc(:,:);
deltasThb(:,Nmu:2*Nmu-1) = deltasThb(:,:);
deltasT  (:,Nmu:2*Nmu-1) = deltasT  (:,:);

deltasc  (:,Nmu-1:-1:1) = conj(deltasc  (:,Nmu+1:2*Nmu-1))*1i;
deltasb  (:,Nmu-1:-1:1) = conj(deltasb  (:,Nmu+1:2*Nmu-1))*1i;
deltasThc(:,Nmu-1:-1:1) = conj(deltasThc(:,Nmu+1:2*Nmu-1))*1i;
deltasThb(:,Nmu-1:-1:1) = conj(deltasThb(:,Nmu+1:2*Nmu-1))*1i;
deltasT  (:,Nmu-1:-1:1) = conj(deltasT  (:,Nmu+1:2*Nmu-1))*1i;

muext              = zeros(1,2*Nmu-1);
muext(Nmu:2*Nmu-1) =  mu(1:Nmu);
muext(Nmu-1:-1:1)  = -mu(2:Nmu);

disp('--- costh between V_cb and wavevector(k) being calculated ---')

%% Calculate norm of V_cb to prevent division by zero
norm_Vcb = norm([V_cb_1_azend(ic,jc,kc) V_cb_2_azend(ic,jc,kc) V_cb_3_azend(ic,jc,kc)]);

if ~memory_save
  if norm_Vcb == 0
    costh_k_V = zeros(size(ksq_p));
  else
    costh_k_V = (V_cb_1_azend(ic,jc,kc)*k1_3D_p + V_cb_2_azend(ic,jc,kc)*k2_3D_p + V_cb_3_azend(ic,jc,kc)*k3_3D_p) / norm_Vcb ./sqrt(ksq_p);
  end
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend = Nmode_p;
    end
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    
    if norm_Vcb == 0
      if kkchunk ~= Nchunk
        costh_k_V(:,:,kkstart:kkend) = zeros(Nmode_p, Nmode_p, Nsub_chunk);
      else
        if Nsub_r == 0
            Nsub = Nsub_chunk;
        else
            Nsub = Nsub_r;
        end
        costh_k_V(:,:,kkstart:kkend) = zeros(Nmode_p, Nmode_p, Nsub);
      end
    else
      if kkchunk ~= Nchunk
        costh_k_V(:,:,kkstart:kkend) = (V_cb_1_azend(ic,jc,kc)*k1_3D_p_chunk ...
              + V_cb_2_azend(ic,jc,kc)*k2_3D_p_chunk ...
              + V_cb_3_azend(ic,jc,kc)*(k3_3D_p_chunk + kshift)) / norm_Vcb ./sqrt(k1k2sq_p_chunk + (k3_3D_p_chunk + kshift).^2);
      else 
        if Nsub_r == 0
          Nsub = Nsub_chunk;
        else
          Nsub = Nsub_r;
        end
        costh_k_V(:,:,kkstart:kkend) = (V_cb_1_azend(ic,jc,kc)*k1_3D_p_chunk(:,:,1:Nsub) ...
              + V_cb_2_azend(ic,jc,kc)*k2_3D_p_chunk(:,:,1:Nsub) ...
              + V_cb_3_azend(ic,jc,kc)*(k3_3D_p_chunk(:,:,1:Nsub) + kshift)) / norm_Vcb ./sqrt(k1k2sq_p_chunk(:,:,1:Nsub) + (k3_3D_p_chunk(:,:,1:Nsub) + kshift).^2);
      end
    end
  end
end

%% Clamp costh_k_V to [-1.0, 1.0] to prevent floating-point errors
%% and nullify the monopole (k=0) NaN.
costh_k_V(costh_k_V > 1.0) = 1.0;
costh_k_V(costh_k_V < -1.0) = -1.0;
costh_k_V(isnan(costh_k_V)) = 0;

%% Generate and dump initial condition data
Generate_IC_CR;  %%==== script ==================

clear costh_k_V;

%% Dump figure-useful data
if matlabflag
  if baryonparticleflag
    save([ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'xbar_plane', 'xbar_ex_plane', 'ybar_plane', 'ybar_ex_plane', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot', '-v6');
  else
    save([ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot', '-v6');
  end
else
  if baryonparticleflag
    save('-mat-binary', [ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'xbar_plane', 'xbar_ex_plane', 'ybar_plane', 'ybar_ex_plane', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot');
  else
    save('-mat-binary', [ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot');
  end
end

disp('*********** bccomics_CR successfully ended ************');
diary off;
movefile('diary', [ICsubdir '/bccomics_CR.log']);
