%% Choose a patch by asking the user for inputs of Delta_c and V_cb.

choose_zi_flag = true; %% By default, choose a patch based on z=1000 condition.

%% get some statistics of patches at z=1000
%% also check isotropy of vector fields which may occur (e.g. velocity)
%% from low-k sampling variance
stdDc   = std(Delta_c(:));
stdTc   = std(Theta_c(:));
stdVc1  = std(V_c_1(:));
stdVc2  = std(V_c_2(:));
stdVc3  = std(V_c_3(:));
rmsVc   = sqrt(mean(V_c_1(:).^2+V_c_2(:).^2+V_c_3(:).^2));
stdDb   = std(Delta_b(:));
stdTb   = std(Theta_b(:));
stdVb1  = std(V_b_1(:));
stdVb2  = std(V_b_2(:));
stdVb3  = std(V_b_3(:));
rmsVb   = sqrt(mean(V_b_1(:).^2+V_b_2(:).^2+V_b_3(:).^2));
stdVcb1 = std(V_cb_1(:));
stdVcb2 = std(V_cb_2(:));
stdVcb3 = std(V_cb_3(:));
rmsVcb  = sqrt(mean(V_cb_1(:).^2+V_cb_2(:).^2+V_cb_3(:).^2));
Vcbp    = sqrt(2/3)*rmsVcb;
stdT    = std(Delta_T(:));
stdcom  = std(Deltacom(:));
stdstr  = std(Deltastr(:));
stdgro  = std(Deltagro(:));
stddec  = std(Deltadec(:));

%% PDF of relV at z=1000
if plotflag
  ifig = ifig+1;
  figure(ifig);
  %% actual data
  if matlabflag
    histogram(Vcb,100,'Normalization','probability');
  else
    hist(Vcb,100,1);
  end
  hold on;
  stdVcb_1D_kms = rmsVcb*MpcMyr_2_kms/sqrt(3)
  relVV = stdVcb_1D_kms*(0:0.05:30);
  %% theoretical
  fVV = sqrt(2/pi)*(relVV.^2/stdVcb_1D_kms^3).*exp(-relVV.^2 /2/stdVcb_1D_kms^2);
  plot(relVV,fVV);
  axis([0 100 0 0.04])
end

%% In case azend values are used...
sDc_azend   = std(Dc3D_azend(:));
Vcb_azend   = sqrt(V_cb_1_azend(:).^2+V_cb_2_azend(:).^2+V_cb_3_azend(:).^2);
rmsVcb_azend=sqrt(mean(Vcb_azend.^2));
Vcbp_azend  = sqrt(2/3)*rmsVcb_azend;  %% peak velocity in Maxwell-Boltzmann distr.

ee          = 1e-2;

if choose_zi_flag
  %% When zi values are used: ------------------------------------------- begin
  disp(['Standard deviation of CDM overdensities (sDc) is ' num2str(stdDc)]);
  disp('Choose CDM overdensity environment:');
  odflag=input('Input 0 for mean, 1 for overdense, 2 for underdense:');
  if (odflag==0)
    odnum=0;
  elseif (odflag==1)
    disp('What multiple of sDc away from the mean overdensity, 0? Example: for Delta_c = +1.5*sDc, Enter 1.5');
    odnum = input('Enter a floating-point number:');
    odnum = abs(odnum)*stdDc; %% into actual value
  elseif (odflag==2)
    disp('What multiple of sDc away from the mean, 0? Example: for Delta_c = -1.5*sDc, Enter 1.5');
    disp('Do not worry about the negative sign, the code knows.');
    odnum = input('Enter a floating-point number: ');
    odnum = -abs(odnum)*stdDc; %% into actual value
  else
    disp('Wrong choice.');
    return;
  end
  disp(['CDM overdensity chosen: Delta_c = ' num2str(odnum) '*sDc = ' num2str(odnum*stdDc)]);
  disp('---------------------------------------');

  disp(['RMS of Vbc (rmsV) at z = ' num2str(zi) ' is ' num2str(rmsVcb*MpcMyr_2_kms) ' km/s']);
  disp(['Peak of Vbc in Maxwell-Boltzmann distribution is ' num2str(Vcbp*MpcMyr_2_kms) ' km/s']);
  disp(['Choose Vbc environment at z = ' num2str(zi)]);
  Vcbnum = input(['Enter Vbc at z = ' num2str(zi) ' in units of km/s: ']);
  Vcbnum = abs(Vcbnum) /MpcMyr_2_kms; %% into Mpc/Myr unit

else
  %% When zzend values are used: ------------------------------------------- begin
  disp(['Standard deviation of CDM overdensities (sDc) is ' num2str(sDc_azend)]);
  disp('Choose CDM overdensity environment:');
  odflag=input('Input 0 for mean, 1 for overdense, 2 for underdense:');
  if (odflag==0)
    odnum=0;
  elseif (odflag==1)
    disp('What multiple of sDc away from the mean overdensity, 0? Example: for Delta_c = +1.5*sDc, Enter 1.5');
    odnum = input('Enter a floating-point number:');
    odnum = abs(odnum)*sDc_azend; %% into actual value
  elseif (odflag==2)
    disp('What multiple of sDc away from the mean, 0? Example: for Delta_c = -1.5*sDc, Enter 1.5');
    disp('Do not worry about the negative sign, the code knows.');
    odnum = input('Enter a floating-point number: ');
    odnum = -abs(odnum)*sDc_azend; %% into actual value
  else
    disp('Wrong choice.');
    return;
  end
  disp(['CDM overdensity chosen: Delta_c = ' num2str(odnum) '*sDc = ' num2str(odnum*sDc_azend)]);
  disp('---------------------------------------');

  disp(['RMS of Vbc (rmsV) at z = ' num2str(zzend) ' is ' num2str(rmsVcb_azend*MpcMyr_2_kms) ' km/s']);
  disp(['Peak of Vbc in Maxwell-Boltzmann distribution is ' num2str(Vcbp_azend*MpcMyr_2_kms) ' km/s']);
  disp(['Choose Vbc environment at z = ' num2str(zzend)]);
  disp(['Note that V_cb is proportional to (1+z)']);
  Vcbnum = input(['Enter Vbc at z = ' num2str(zzend) ' in units of km/s: ']);
  Vcbnum = abs(Vcbnum) /MpcMyr_2_kms; %% into Mpc/Myr unit

  %% Find index of patches with user-selected overdensity and Vcb (with ~1% margin)
  %% Need ee*sDc and ee*rmsVcb for odnum=0 case.
  %% numeric flags are multiplied below to mimic "AND" boolean
  ind_od  = find(((1-ee)*odnum -ee*sDc_azend <= Dc3D_azend(:)).*(Dc3D_azend(:) <= (1+ee)*odnum +ee*sDc_azend   ));
  ind_vcb = find(((1-ee)*Vcbnum              <= Vcb_azend(:) ).*(Vcb_azend(:)  <= (1+ee)*Vcbnum+ee*rmsVcb_azend));

  %% indices of patches satisfying both conditions
  indices_patch = ind_od(ismember(ind_od,ind_vcb));
  if (length(indices_patch)>0)
    disp([num2str(length(indices_patch)) ' patches out of total ' num2str(Ncell^3) ' patches satisfy your chosen condition with 1% margin.']);
  else %% In case no patch is found, relax the condition
    disp('Loosening patch finding condition to 2%');
    ind_od  = find(((1-2*ee)*odnum -2*ee*sDc_azend <= Dc3D_azend(:)).*(Dc3D_azend(:) <= (1+2*ee)*odnum +2*ee*sDc_azend   ));
    ind_vcb = find(((1-2*ee)*Vcbnum                <= Vcb_azend(:) ).*(Vcb_azend(:)  <= (1+2*ee)*Vcbnum+2*ee*rmsVcb_azend));

    indices_patch = ind_od(ismember(ind_od,ind_vcb));
    if (length(indices_patch)==0)
      disp('No such patch exists. Note that Vcb=0 case is very rare by nature!!');
      disp('Also check whether Vcb is in km/s and not like 6 sigma away from zero.');
      return;  
    end
  end

  %% Find the best matching patch
  disp('----------------One best matching patch is being found----------------');
  relerr = (Dc3D_azend(indices_patch)-odnum).^2/sDc_azend^2 + (Vcb_azend(indices_patch)-Vcbnum).^2/rmsVcb_azend^2;
  [minrelerr, indmin] = min(relerr);
  ind_patch = indices_patch(indmin);
  disp(['Wanted Delta_c = ' num2str(odnum*sDc_azend) '; Selected patch''s Delta_c = ' num2str(Dc3D_azend(ind_patch))]);
  disp(['Wanted Vcb = ' num2str(Vcbnum*MpcMyr_2_kms) ' km/s; Selected patch''s Vcb = ' num2str(Vcb_azend(ind_patch)*MpcMyr_2_kms) ' km/s']);
  %% When azend values are used: --------------------------------------------- end
end

%% Convert 1D index into 3D index
[icc1 icc2 icc3] = ind2sub([Nmode, Nmode, Nmode], ind_patch);
icc = [icc1 icc2 icc3];

%% Check if your patch has been selected already
if exist([outputdir '/icc.dat'], 'file')
  icc_read=load([outputdir '/icc.dat']);
  if ismember(icc,icc_read)
    disp('The patch has already been selected. Quitting');
    return;
  end
end

fout=fopen([outputdir '/icc.dat'],'a');
fprintf(fout,'%i %i %i\n',icc');
fclose(fout);

daticc(1) = Dc3D_azend  (icc1,icc2,icc3);
daticc(2) = Db3D_azend  (icc1,icc2,icc3);
daticc(3) = THc3D_azend (icc1,icc2,icc3);
daticc(4) = THb3D_azend (icc1,icc2,icc3);
daticc(5) = V_cb_1_azend(icc1,icc2,icc3)*MpcMyr_2_kms;
daticc(6) = V_cb_2_azend(icc1,icc2,icc3)*MpcMyr_2_kms;
daticc(7) = V_cb_3_azend(icc1,icc2,icc3)*MpcMyr_2_kms;
daticc(8) = sqrt(daticc(5).^2 + daticc(6).^2 + daticc(7).^2);
daticc(9) = DT3D_azend  (icc1,icc2,icc3);
fout=fopen([outputdir '/icc_Dc_Db_Thc_Thb_Vcb1_Vcb2_Vcb3_Vcb_DT.dat'],'a');
fprintf(fout,'%i %i %i %e %e %e %e %e %e %e %e %e\n',[icc daticc]');
fclose(fout);

