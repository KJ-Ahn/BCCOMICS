%%
%% BCCOMICS: Reads in fluctuation (curvature*TF) made by bccomics_setup.m.
%%           Generates 3D fields of small-scale (inside-a-patch)
%%           perturbations into bare binary files. 
%%
%% Author: Kyungjin Ahn
%%
%% This MATLAB(R) / GNU Octave code is freely distributed, and you are
%% free to modify it or port it into other languages. BCCOMICS is under
%% an absolutely no-warranty condition. It is assumed that you consent
%% to one condition: when you get scientific results using BCCOMICS
%% and publish them, in your paper you need to cite this paper
%% (please use journal-provided id after publication),
%% ---------
%% Ahn & Smith 2018, arXiv:1807.04063 (AS18).
%% ---------
%% For detailed theoretical background, please cite this paper
%% (not a requirement for using this code though)
%% ---------
%% Ahn 2016, ApJ 830:68 (A16).
%% ---------
%%
%% Other references:
%%   Ma & Bertschinger 1995, ApJ 455, 7 (MB)
%%   Naoz & Barkana 2005, MNRAS 362, 1047 (NB)
%%   Tseliakhovich & Hirata 2010, PRD, 82, 083520 (TH)
%%
%%
%% What it does: This code reads in a fluctuation (~curvature*transfer-function)
%%               and monopole (e.g. V_cb of a patch) data generated by 
%%               bccomics_setup. For any wavevector k, the read in data
%%               is interpolated onto (k, mu). Then random seed is applied,
%%               and FFT is performed.
%%               For initial condition readable by enzo, a conversion script
%%               of this binary into hdf5 is provided ("convert_enzo.py").
%%               Conversion sripts for other codes are welcomed!!
%% 
%%               
%% Some details:
%% ----------
%% Just uniform-grid initial condition only. No nested grid IC yet.
%% Binary files from bccomics_setup.m is also in hdf5 format, so
%% porting bccomics.m into other languages & improving it are welcomed.

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
%% index for center of k-space (e.g. if 6 sample points exist, 4th is the
%% center, not 3rd). This convention for even number is different from that
%% in p.69 of "DFT: An Owner Manual ..." by W. Briggs.
%% k index runs from -N/2 to N/2-1 in this code, but Briggs uses
%% -N/2+1 to N/2. Had to choose the former convention due to FFT convention
%% of Matlab and Octave for even numbered cases.

%% Read in parameters for initial condition
patch_init;  %%==== script ==================

interp2opt = 'cubic'
%% May choose 'pchip' for Matlab below, but for consistency with Octave
%% just use 'linear'. Octave interpn does not have 'pchip' implemented yet.
%% 'spline' is somewhat dangerous.
interpnopt = 'linear' 

%% For assigning k, see p.69 of "DFT..." by W. Briggs.
%% The convention below has [-Nhalf_p:Nhalf_p-1], different
%% from Briggs convention [-Nhalf_p+1:Nhalf_p], but this is
%% to par with Matlab and Octave FFT convention.


%% k1 component on each (k1,k2,k3) point, etc.
%% If memory error occurs, better to turn on memory_save (looping over k3 axis)
memory_save = true;
if Nmode_p<=32
  memory_save = false;
end

%% number of chunks to iterate for memory_save
if memory_save
  Nsub_chunk = 8;  %% # of slices to treat in one chunk
  Nchunk = ceil(Nmode_p/Nsub_chunk);  %% check for residual slices after chunkening
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
  %% makes Nmode_p*Nmode_p*Nsub_chunk k1-vector & k2-vector array, using cat(concatenate).
  for kkk=2:Nsub_chunk
    k1_3D_p_chunk = cat(3, k1_3D_p_chunk, k1_2D_p);
    k2_3D_p_chunk = cat(3, k2_3D_p_chunk, k2_2D_p);
  end

  onesk1k2_2D = ones(Nmode_p,Nmode_p);  %% ones on k1-k2 plane
  k3_3D_p_chunk = kunit_p*(-Nhalf_p)*onesk1k2_2D;  %% 2D array at this point
  %% makes Nmode_p*Nmode_p*Nsub_chunk k1-vector & k2-vector array, using cat(concatenate).
  %% integer index runs from -Nhalf_p:-Nhalf_p+(Nsub_chunk-1)
  for kkk=-Nhalf_p+1:-Nhalf_p+(Nsub_chunk-1)
    k3_3D_p_chunk = cat(3, k3_3D_p_chunk, kunit_p*kkk*onesk1k2_2D);
  end
  
  k1k2sq_p_chunk = k1_3D_p_chunk.^2 +k2_3D_p_chunk.^2; %% k1^2+k2^2
end

%% utilize above for rvector too, but just in memory saving way (****)
%%r1 = k1_3D_p/kunit_p;
%%r2 = k2_3D_p/kunit_p;
%%r3 = k3_3D_p/kunit_p;

%% read in mu info
mu  = load([setupdir '/mu.dat']);
dmu = mu(2)-mu(1);
Nmu = length(mu);

%% choose patch to generate initial condition on
Choose_finalpatch;  %%==== script ==================
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

ICsubdir = [ICdir '/' num2str(Lbox_p_inMpch,'%.2f') 'Mpch_' num2str(Ncell_p) '_ic' num2str(ic) '_jc' num2str(jc) '_kc' num2str(kc)];
if ~exist(ICsubdir)
  mkdir(ICsubdir); 
end

%% Copy essential files to ICsubdir. This helps consistency in record tracking.
copyfile('params.m', ICsubdir);
copyfile('params_patch.m', ICsubdir);
copyfile(Cosmology, ICsubdir);
copyfile([setupdir '/zz.dat'], ICsubdir);
copyfile([setupdir '/stats_zend.dat'], ICsubdir);
iccdat = cellspec_azend(idxcc,:);
fout   = fopen([ICsubdir '/icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat'],'w');
fprintf(fout, '%4i %4i %4i %e %e %e %e %e %e %e %e %e\n', iccdat');
fclose(fout);

zf = zzend;  %% redshift for initial condition
af = 1/(1+zf);  %% scale factor for initial condition

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
  %% check if wrongful (in size) matbin file is written, and if so write in simpler binary.
  %% (Octave behaves badly when writing matbin file with big (~1000*1000*500*2)array)
  if matlabflag
    msg   = dir(fileNseed);
    fsize = msg.bytes; %% file size in bytes
  else
    msg   = lstat(fileNseed);
    fsize = msg.size; %% file size in bytes
  end  
  if (fsize < Nmode_p*Nmode_p*Nc_p*2*8)  %% randamp & randphs at 8 bytes
    delete(fileNseed);
    fileNseed_1 = [ICsubdir '/subgaussseed' num2str(Nmode_p) '.bin'];
    ffout = fopen(fileNseed_1,'w');
    fwrite(ffout, randamp, 'double');
    fwrite(ffout, randphs, 'double');
    fclose(ffout);
  end
end


%% For a given k, -mu case has its Real same as Imag of mu case,
%%                         and its Imag same as Real of mu case.
%% Switching Real and Imag is done easily by i*conj(complex_number).
%% -- First, shift mu=[0,...,1] values to right (matrices increase in size).
deltasc  (:,Nmu:2*Nmu-1) = deltasc  (:,:);
deltasb  (:,Nmu:2*Nmu-1) = deltasb  (:,:);
deltasThc(:,Nmu:2*Nmu-1) = deltasThc(:,:);
deltasThb(:,Nmu:2*Nmu-1) = deltasThb(:,:);
deltasT  (:,Nmu:2*Nmu-1) = deltasT  (:,:);
%% -- Then, generate mu=[-1,...,0) values
deltasc  (:,Nmu-1:-1:1) = conj(deltasc  (:,Nmu+1:2*Nmu-1))*i;
deltasb  (:,Nmu-1:-1:1) = conj(deltasb  (:,Nmu+1:2*Nmu-1))*i;
deltasThc(:,Nmu-1:-1:1) = conj(deltasThc(:,Nmu+1:2*Nmu-1))*i;
deltasThb(:,Nmu-1:-1:1) = conj(deltasThb(:,Nmu+1:2*Nmu-1))*i;
deltasT  (:,Nmu-1:-1:1) = conj(deltasT  (:,Nmu+1:2*Nmu-1))*i;
  
%% Extend mu to cover full angle accordingly: muext=[-1,...,0,...,1]
muext              = zeros(1,2*Nmu-1);
muext(Nmu:2*Nmu-1) =  mu(1:Nmu);
muext(Nmu-1:-1:1)  = -mu(2:Nmu);

%% mu = cosine(angle between k vector and V_cb=V_c-V_b).
%% The mu convention is consistent with f.m.
disp('--- costh between V_cb and wavevector(k) being calculated ---')
if ~memory_save
  costh_k_V = (V_cb_1_azend(ic,jc,kc)*k1_3D_p + V_cb_2_azend(ic,jc,kc)*k2_3D_p + V_cb_3_azend(ic,jc,kc)*k3_3D_p) /norm([V_cb_1_azend(ic,jc,kc) V_cb_2_azend(ic,jc,kc) V_cb_3_azend(ic,jc,kc)]) ./sqrt(ksq_p);
else
  for kkchunk=1:Nchunk
    disp(['*** ' num2str(kkchunk) ' out of ' num2str(Nchunk) ' chunks being processed ***'])
    kkstart = 1 + (kkchunk-1)*Nsub_chunk;
    kkend   = kkchunk*Nsub_chunk;
    if kkchunk == Nchunk
      kkend = Nmode_p;
    end
    %% k3_3D_p_chunk is defined at the bottom chunk, so need to add a shift in iteration
    kshift = kunit_p*(kkchunk-1)*Nsub_chunk;
    if kkchunk ~= Nchunk
      costh_k_V(:,:,kkstart:kkend) = (V_cb_1_azend(ic,jc,kc)*k1_3D_p_chunk             + V_cb_2_azend(ic,jc,kc)*k2_3D_p_chunk             + V_cb_3_azend(ic,jc,kc)*(k3_3D_p_chunk            +kshift)) /norm([V_cb_1_azend(ic,jc,kc) V_cb_2_azend(ic,jc,kc) V_cb_3_azend(ic,jc,kc)]) ./sqrt(k1k2sq_p_chunk            +(k3_3D_p_chunk            +kshift).^2);
    else %% the last chunk
      if Nsub_r == 0
        Nsub = Nsub_chunk;
      else
        Nsub = Nsub_r;
      end
      costh_k_V(:,:,kkstart:kkend) = (V_cb_1_azend(ic,jc,kc)*k1_3D_p_chunk(:,:,1:Nsub) + V_cb_2_azend(ic,jc,kc)*k2_3D_p_chunk(:,:,1:Nsub) + V_cb_3_azend(ic,jc,kc)*(k3_3D_p_chunk(:,:,1:Nsub)+kshift)) /norm([V_cb_1_azend(ic,jc,kc) V_cb_2_azend(ic,jc,kc) V_cb_3_azend(ic,jc,kc)]) ./sqrt(k1k2sq_p_chunk(:,:,1:Nsub)+(k3_3D_p_chunk(:,:,1:Nsub)+kshift).^2);
    end
  end
end

%% Generate and dump initial condition data
Generate_IC;  %%==== script ==================

clear costh_k_V;

%% Dump figure-useful data
if matlabflag
  if baryonparticleflag
    save([ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'xbar_plane', 'xbar_ex_plane', 'ybar_plane', 'ybar_ex_plane', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot', '-v6')
  else
    save([ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot', '-v6')
  end
else
  if baryonparticleflag
    save('-mat-binary', [ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'xbar_plane', 'xbar_ex_plane', 'ybar_plane', 'ybar_ex_plane', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot')
  else
    save('-mat-binary', [ICsubdir '/4fig.matbin'], 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot')
  end
end

disp('*********** bccomics successfully ended ************');
diary off;
movefile('diary', [ICsubdir '/bccomics.log']);

