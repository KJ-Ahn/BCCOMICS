%% This script generates initial condition data for enzo and dumps them.
%% For ICs for other simulation codes, this is the script to start from.
%%
%% Initial condition generation using following steps, for enzo.
%% (1) Get k and mu dependent fluctuation using transfer function.
%% (2) Apply random seed and normalize (rand_real_norm).
%% (3) FFT & record.


%% =========== CDM density and position ======================== begin

%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
%% with array size Nmode_p, Nmode_p, Nmode_p.
%% Matlab allows extrapolation only for 'spline' method.
%% (costh_k_V and ksq_p have same array dimension, so interp2 
%%  interprets these as scattered data points: see interp2 instruction)
%% This is linear logarithmic interpolation along k, so the monopole
%% term (k=0) may obtain inf or nan due to 0.5*log(ksq_p). 
%% We will cure this by nullifying monopole anyway down below (**).
dc  = interp2(muext,log(ksampletab), deltasc,  costh_k_V,0.5*log(ksq_p),interp2opt);  %% dc still k-space values here.

%% randomize, apply reality, and normalize
dc = rand_real_norm(dc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);
%% CDM displacement vector, related to CDM density at 1st order.
%% No need for above normalization because this is
%% derived after above normalization on dc.
%% ------------- cpos1 ----------------------
Psi1                 = i*k1_3D_p./ksq_p.*dc;
Psi1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi1                 = real(ifftn(ifftshift(Psi1)));
%% Normalized position of particles in domain [0,1), cell-centered way. (enzo)
%% If unperturbed(Psi=0), it should run [0.5, 1.5, ...., Nmode_p-0.5]/Nmode_p,
%% For enzo, wrapping needed if perturbed potition is out of the domain [0, 1).
fout = fopen([ICsubdir '/cpos1'], 'w');
fwrite(fout, mod((Psi1 + (k1_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
fclose(fout);
xCDM_plane    =   Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure
xCDM_ex_plane = 5*Psi1(:,:,1) + k1_3D_p(:,:,1)/kunit_p*Lcell_p + Lbox_p/2; %% for figure, NOT REAL but to make more contrast in CDM position
%% Prepare for interpn, let positions run from 1:Nmode_p for unperturbed particles,
%% to get more-accurate-than-1LPT velocity when particlevelocity_accuracyflag=true.
%% Using mod function, the actual positions will run from 1 to Nmode_p+0.9999999...
%% Interpolation basis will have domain 1:Nmode_p+1 for safe interpolation (see **).
if particlevelocity_accuracyflag  
  Psi1 = mod((Psi1 + (k1_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
else
  clear Psi1  %% save memory
end

%% ------------- cpos2 ----------------------
Psi2                 = i*k2_3D_p./ksq_p.*dc;
Psi2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi2                 = real(ifftn(ifftshift(Psi2)));
fout = fopen([ICsubdir '/cpos2'], 'w');
fwrite(fout, mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
fclose(fout);
if particlevelocity_accuracyflag
  Psi2 = mod((Psi2 + (k2_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
else
  clear Psi2  %% save memory
end

%% ------------- cpos3 ----------------------
Psi3                 = i*k3_3D_p./ksq_p.*dc;
Psi3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
Psi3                 = real(ifftn(ifftshift(Psi3)));
fout = fopen([ICsubdir '/cpos3'], 'w');
fwrite(fout, mod((Psi3 + (k3_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
fclose(fout);
if particlevelocity_accuracyflag
  Psi3 = mod((Psi3 + (k3_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
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
Thc = interp2(muext,log(ksampletab), deltasThc, costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
Thc = rand_real_norm(Thc,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vc1 ----------------------
vc1(:,:,:)          = -i*af*k1_3D_p./ksq_p.*Thc(:,:,:);
vc1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc1                 = real(ifftn(ifftshift(vc1)));
if particlevelocity_accuracyflag
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc1 = padarray(vc1, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc1 = interpn(vc1, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
%% enzo velocity output is the following:
%%vc1_enzo = vc1 * MpcMyr_2_kms * 1e5 /VelocityUnits;
fout = fopen([ICsubdir '/vc1'], 'w');
fwrite(fout, vc1*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);
Vc1 = reshape(vc1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vc1  %% save memory

%% ------------- vc2 ----------------------
vc2(:,:,:)          = -i*af*k2_3D_p./ksq_p.*Thc(:,:,:);
vc2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc2                 = real(ifftn(ifftshift(vc2)));
if particlevelocity_accuracyflag
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc2 = padarray(vc2, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc2 = interpn(vc2, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
fout = fopen([ICsubdir '/vc2'], 'w');
fwrite(fout, vc2*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);
Vc2 = reshape(vc2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vc2  %% save memory

%% ------------- vc3 ----------------------
vc3(:,:,:)          = -i*af*k3_3D_p./ksq_p.*Thc(:,:,:);
vc3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vc3                 = real(ifftn(ifftshift(vc3)));
if particlevelocity_accuracyflag
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vc3 = padarray(vc3, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vc3 = interpn(vc3, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
fout = fopen([ICsubdir '/vc3'], 'w');
fwrite(fout, vc3*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);
Vc3 = reshape(vc3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vc3  %% save memory

%% ------------- velocity divergence ---------------------
Thc = real(ifftn(ifftshift(Thc)));  

ZThc = reshape(Thc(:,:,1),Nmode_p,Nmode_p); %% for figure
clear Thc  %% save memory
%% =========== CDM velocity ==================================== end



%% =========== baryon density ================================== begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
db  = interp2(muext,log(ksampletab), deltasb,  costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
db = rand_real_norm(db,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%%%% If SPH particle is used, one can here get the particle positions 
%%%% just the way CDM positions are calculated.
if baryonparticleflag
  %% ------------- bpos1 ----------------------
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
  Psi2                 = i*k2_3D_p./ksq_p.*db;
  Psi2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
  Psi2                 = real(ifftn(ifftshift(Psi2)));
  fout = fopen([ICsubdir '/bpos2'], 'w');
  fwrite(fout, mod((Psi2 + (k2_3D_p/kunit_p+0.5)*Lcell_p + Lbox_p/2)/Lbox_p, 1), 'double');
  fclose(fout);
  if particlevelocity_accuracyflag
    Psi2 = mod((Psi2 + (k2_3D_p/kunit_p)*Lcell_p + Lbox_p/2)/Lbox_p*Nmode_p, Nmode_p)+1;
  else
    clear Psi2  %% save memory
  end
  
  %% ------------- bpos3 ----------------------
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
db = real(ifftn(ifftshift(db)));  

%% enzo baryon density output is the following:
%%  db_enzo    = (db+1)*fb
fout = fopen([ICsubdir '/db'], 'w');
fwrite(fout, (db+1)*fb, 'double');
fclose(fout);

Zb    = reshape(db(:,:,1),Nmode_p,Nmode_p);
clear db  %% save memory
%% =========== baryon density ================================== end



%% =========== baryon velocity ==================================== begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
Thb = interp2(muext,log(ksampletab), deltasThb, costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
Thb = rand_real_norm(Thb,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- vb1 ----------------------
vb1(:,:,:)          = -i*af*k1_3D_p./ksq_p.*Thb(:,:,:);
vb1(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb1                 = real(ifftn(ifftshift(vb1)));
%% add streaming velocity (V_cb = Vc - Vb)
vb1                 = vb1 - V_cb_1_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--1--**)
sp_Etot_enzo = 1/2*(vb1*MpcMyr_2_kms*1e5/VelocityUnits).^2; 
if (particlevelocity_accuracyflag & baryonparticleflag)
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb1 = padarray(vb1, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb1 = interpn(vb1, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
%% enzo velocity output is the following:
%%vb1_enzo = vb1 * MpcMyr_2_kms * 1e5 /VelocityUnits;
fout = fopen([ICsubdir '/vb1'], 'w');
fwrite(fout, vb1*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);

Vb1 = reshape(vb1(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vb1  %% save memory

%% ------------- vb2 ----------------------
vb2(:,:,:)          = -i*af*k2_3D_p./ksq_p.*Thb(:,:,:);
vb2(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb2                 = real(ifftn(ifftshift(vb2)));
%% add streaming velocity (V_cb = Vc - Vb)
vb2                 = vb2 - V_cb_2_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--2--**)
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb2*MpcMyr_2_kms*1e5/VelocityUnits).^2; 
if (particlevelocity_accuracyflag & baryonparticleflag)
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb2 = padarray(vb2, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb2 = interpn(vb2, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
fout = fopen([ICsubdir '/vb2'], 'w');
fwrite(fout, vb2*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);

Vb2 = reshape(vb2(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vb2  %% save memory

%% ------------- vb3 ----------------------
vb3(:,:,:)          = -i*af*k3_3D_p./ksq_p.*Thb(:,:,:);
vb3(Nc_p,Nc_p,Nc_p) = complex(0);  %% fixing nan or inf monopole
vb3                 = real(ifftn(ifftshift(vb3)));
%% add streaming velocity (V_cb = Vc - Vb)
vb3                 = vb3 - V_cb_3_azend(ic,jc,kc); 
%% memory-saving way of calculating sp_Etot_enzo (**--3--**)
sp_Etot_enzo = sp_Etot_enzo + 1/2*(vb3*MpcMyr_2_kms*1e5/VelocityUnits).^2; 
if (particlevelocity_accuracyflag & baryonparticleflag)
  %% pad with periodic boundary condition, to provide 1:Nmode_p+1 domain. (**)
  vb3 = padarray(vb3, [1 1 1], 'circular', 'post'); %% now (Nmode_p+1)^3 elements.
  %% Find 
  vb3 = interpn(vb3, Psi1, Psi2, Psi3, interpnopt); %% now Nmode_p^3 elements
end
fout = fopen([ICsubdir '/vb3'], 'w');
fwrite(fout, vb3*MpcMyr_2_kms*1e5/VelocityUnits, 'double');
fclose(fout);

Vb3 = reshape(vb3(:,:,1) *MpcMyr_2_kms, Nmode_p, Nmode_p); %% for figure
clear vb3  %% save memory

%% ------------- velocity divergence ---------------------
Thb = real(ifftn(ifftshift(Thb)));  

ZThb = reshape(Thb(:,:,1),Nmode_p,Nmode_p); %% for figure
clear Thb  %% save memory
%% =========== baryon velocity ==================================== end


%% =========== baryon temperature, energies ======================= begin
%% Matlab & Octave 2D interpolation!! --> generating k-space deltas 
dT  = interp2(muext,log(ksampletab), deltasT,  costh_k_V,0.5*log(ksq_p),interp2opt);

%% randomize, apply reality, and normalize
dT = rand_real_norm(dT,Nmode_p,Nc_p,randamp,randphs,Vbox_p);

%% ------------- temperature, energies  ---------------------
dT = real(ifftn(ifftshift(dT)));  


%% Mean IGM temperature fit from Tseliakhovich & Hirata
aa1  = 1/119
aa2  = 1/115
Tz    = TCMB0/af /(1+af/aa1/(1+(aa2/af)^1.5));  %% in K, global average temperature
Tz    = Tz*(1+DT3D_azend(ic,jc,kc)); %% local(cell) average temperature
  
%% specific thermal energy for monatomic gas (H+He), in units of VelocityUnits^2 : sp_Eth_enzo
%%Tcell = (dT+1)*Tz;  %% in K
%%sp_Eth_enzo  = 3/2*kb*Tcell /(mmw*mH) /VelocityUnits^2;  %% see Enzo paper(2014) eq. 7.
fout = fopen([ICsubdir '/etherm'], 'w');
fwrite(fout, 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2, 'double');
fclose(fout);

%% Zeth in erg (thermal energy per baryon)
Zeth = reshape(3/2*kb*(dT(:,:,1)+1)*Tz/(mmw*mH),Nmode_p,Nmode_p); %% for figure

%% Ztemp in K
Ztemp = reshape((dT(:,:,1)+1)*Tz,Nmode_p,Nmode_p); %% for figure


%% specific total energy for monatomic gas (H+He), in units of VelocityUnits^2 :   sp_Etot_enzo
%% Currently sp_Etot_enzo does not include magnetic contribution, but in principle it should.
%% memory-saving way of calculating sp_Etot_enzo (**--4--**)
sp_Etot_enzo = sp_Etot_enzo + 3/2*kb*(dT+1)*Tz /(mmw*mH) /VelocityUnits^2;
  
fout = fopen([ICsubdir '/etot'], 'w');
fwrite(fout, sp_Etot_enzo, 'double');
fclose(fout);

%% Zetot in erg (total energy per baryon)
Zetot = reshape(sp_Etot_enzo(:,:,1)*VelocityUnits^2,Nmode_p,Nmode_p); %% for figure

clear dT sp_Etot_enzo
%% =========== baryon temperature, energies ======================= begin

%% Save some memory
clear costh_k_V 
%clear randamp randphs k1_3D_p k2_3D_p k3_3D_p ksq_p
