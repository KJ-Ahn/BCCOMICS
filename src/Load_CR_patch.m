%% Script: Load_CR_patch.m
%% What it does: Replaces Choose_finalpatch.m for the CR pipeline.
%%               Automatically loads the CR target patch from accumulated records
%%               and prompts the user for selection with the original UI style.

if matlabflag
  load([setupdir '/V_cb_1_azend.matbin'], '-mat', 'V_cb_1_azend');
  load([setupdir '/V_cb_2_azend.matbin'], '-mat', 'V_cb_2_azend');
  load([setupdir '/V_cb_3_azend.matbin'], '-mat', 'V_cb_3_azend');
  load([setupdir '/DT_azend.matbin'],     '-mat', 'DT3D_azend');
  load([setupdir '/Dc3D_azend.matbin'],   '-mat', 'Dc3D_azend');
  load([setupdir '/Db3D_azend.matbin'],   '-mat', 'Db3D_azend');
  load([setupdir '/THc3D_azend.matbin'],  '-mat', 'THc3D_azend');
  load([setupdir '/THb3D_azend.matbin'],  '-mat', 'THb3D_azend');
else
  load('-mat-binary', [setupdir '/V_cb_1_azend.matbin'], 'V_cb_1_azend');
  load('-mat-binary', [setupdir '/V_cb_2_azend.matbin'], 'V_cb_2_azend');
  load('-mat-binary', [setupdir '/V_cb_3_azend.matbin'], 'V_cb_3_azend');
  load('-mat-binary', [setupdir '/DT_azend.matbin'],     'DT3D_azend');
  load('-mat-binary', [setupdir '/Dc3D_azend.matbin'],   'Dc3D_azend');
  load('-mat-binary', [setupdir '/Db3D_azend.matbin'],   'Db3D_azend');
  load('-mat-binary', [setupdir '/THc3D_azend.matbin'],  'THc3D_azend');
  load('-mat-binary', [setupdir '/THb3D_azend.matbin'],  'THb3D_azend');
end

cellspec = load([setupdir '/zi_icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat']);
cellspec_azend = load([setupdir '/zend_icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat']);
Ncc = size(cellspec, 1);

%% Read in z=zi=1000 statistics to normalize Delta_m
fin=fopen([setupdir '/stats_zi.dat']);
fgets(fin); %% skip a line
fgets(fin); %% skip another line
statszi = fscanf(fin, '%e %e %e %e %e %e %e %e %e %e');
fclose(fin);

%% Let user choose a patch with the original UI style + Relative Direction
disp('Patches ordered in calculation time, from oldest(top) to newest(bottom)');
disp('---------------------------------------------------------------------------------------');
disp('Patch #  ix  iy  iz  Delta_m/sigma(Delta_m)   V_bc(km/s)    Vbc_dir[x,y,z]   at z=1000');
for ip=1:Ncc
  Dm_ip = fc * cellspec(ip,4) + fb * cellspec(ip,5);
  Vmag  = cellspec(ip,11);
  
  %% cellspec(8:10) stores V_cb. To display V_bc properly, we invert the signs.
  Vdir_raw = -cellspec(ip,8:10); 
  
  if Vmag > 1e-6
      %% Find the minimum non-zero component to restore relative integer ratios
      nonzero_V = abs(Vdir_raw(abs(Vdir_raw) > 1e-6));
      Vdir_rel = Vdir_raw / min(nonzero_V);
      Vdir_rel = round(Vdir_rel * 10) / 10; %% Clean up floating-point errors
  else
      Vdir_rel = [0, 0, 0];
  end
  AA = [ip cellspec(ip,1) cellspec(ip,2) cellspec(ip,3) Dm_ip/statszi(10) Vmag Vdir_rel(1) Vdir_rel(2) Vdir_rel(3)];
  fprintf('%3i     %3i %3i %3i     %10.3e           %10.3e     [%g, %g, %g]\n',AA);
end

disp(['Choose a patch of your interest; default is ' num2str(Ncc) ' if you just hit Enter below.']);
idxcc = input('Enter your choice (patch #):');
if isempty(idxcc)  %% default to the last patch calculated
  idxcc=Ncc;
end
disp(['Patch # ' num2str(idxcc) ' chosen.']);

ic   = cellspec(idxcc,1);
jc   = cellspec(idxcc,2);
kc   = cellspec(idxcc,3);

%% --- SAFETY GUARD: Prevent loading overwritten patches ---
r   = cellspec_azend(idxcc,:);
cur = [Dc3D_azend(ic,jc,kc) Db3D_azend(ic,jc,kc) ...
       [V_cb_1_azend(ic,jc,kc) V_cb_2_azend(ic,jc,kc) V_cb_3_azend(ic,jc,kc)]*MpcMyr_2_kms];
rec = r([4 5 8 9 10]);
if any(abs(cur - rec) > 1e-5*abs(rec) + 1e-10)
  disp('*** ERROR: This row was overwritten by a later bccomics_CR_setup run.');
  disp('*** Re-run bccomics_CR_setup with this target to generate fresh 3D fields.');
  returnflag = true; return;
end

%% Flush output buffer for Octave to prevent hanging display
if ~matlabflag, fflush(stdout); end
