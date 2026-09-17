%% Script: Set_constraints.m
%% What it does: Replaces Choose_patch.m. Fixes target Delta_m and V_cb
%%               at the center of the box for Constrained Realization.

disp('-----------------------------------------------------------------------');
disp(' Constrained Realization (CR) Setup');
disp('-----------------------------------------------------------------------');

%% Target Overdensity
disp(['Standard deviation of Matter overdensities (stdDm) is ' num2str(stdDm)]);
disp('Choose Matter overdensity environment: ');
odflag=input('Input 0 for mean, 1 for overdense, 2 for underdense: ');

if (odflag==0)
  target_Dm=0;
  odnum_print=0;
elseif (odflag==1)
  disp('What multiple of stdDm away from the mean overdensity, 0? Example: for Delta_m = +1.5*stdDm, Enter 1.5');
  odnum = input('Enter a floating-point number: ');
  target_Dm = abs(odnum)*stdDm; 
  odnum_print = abs(odnum);
elseif (odflag==2)
  disp('What multiple of stdDm away from the mean, 0? Example: for Delta_m = -1.5*stdDm, Enter 1.5');
  disp('Do not worry about the negative sign, the code knows.');
  odnum = input('Enter a floating-point number: ');
  target_Dm = -abs(odnum)*stdDm; 
  odnum_print = -abs(odnum);
else
  disp('Wrong choice.');
  returnflag=true;
  return;
end

disp(['Matter overdensity chosen: Delta_m = ' num2str(odnum_print) '*stdDm = ' num2str(target_Dm)]);
disp('---------------------------------------');

%% Target Velocity Magnitude
Vcbp = sqrt(2/3)*rmsVcb; 
disp(['RMS of Vbc (rmsV) at z = ' num2str(zi) ' is ' num2str(rmsVcb*MpcMyr_2_kms) ' km/s']);
disp(['Peak of Vbc in Maxwell-Boltzmann distribution is ' num2str(Vcbp*MpcMyr_2_kms) ' km/s']);
disp(['Choose Vbc environment at z = ' num2str(zi)]);
vcb_mag_input  = input(['Enter Vbc magnitude at z = ' num2str(zi) ' in units of km/s: ']);
target_Vcb_mag = abs(vcb_mag_input) / MpcMyr_2_kms; 
disp('---------------------------------------');

%% Target Velocity Direction
disp('Enter the velocity direction vector [x, y, z] as an array.');
disp('Example: [1, 1, 1] for isotropic, [1, 0, 0] for x-axis only, or [1, -1, 1].');
v_dir = input('Enter an array: ');

if norm(v_dir) == 0
    v_dir_norm = [0, 0, 0];
else
    v_dir_norm = v_dir / norm(v_dir);
end

target_Vcb_1 = target_Vcb_mag * v_dir_norm(1);
target_Vcb_2 = target_Vcb_mag * v_dir_norm(2);
target_Vcb_3 = target_Vcb_mag * v_dir_norm(3);

%% Set target coordinates (Center of the box)
icc1 = floor(Nmode/2) + 1;
icc2 = floor(Nmode/2) + 1;
icc3 = floor(Nmode/2) + 1;
icc  = [icc1, icc2, icc3];

%% Construct the constraint vector C
%% Note: V_cb constraint is removed from CR filter to prevent k>0 ringing.
%% V_cb will be forced as a background bulk flow (k=0) in Apply_CR_filter.m
C_target = [target_Dm];

disp('----------------Target parameters fixed----------------');
disp(['Target coordinates fixed at center: [' num2str(icc) ']']);
disp(['Target Delta_m: ' num2str(target_Dm)]);
disp(['Target V_bc [x, y, z] in km/s: [' num2str(target_Vcb_1*MpcMyr_2_kms) '  ' num2str(target_Vcb_2*MpcMyr_2_kms) '  ' num2str(target_Vcb_3*MpcMyr_2_kms) ']']);
disp('-------------------------------------------------------');
