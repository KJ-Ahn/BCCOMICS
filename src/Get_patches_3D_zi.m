%% Generate 3D real-space maps at z=zi=1000

k1_3D = zeros(Nmode,Nmode,Nmode);
k2_3D = zeros(Nmode,Nmode,Nmode);
k3_3D = zeros(Nmode,Nmode,Nmode);

%% k1 component on each (k1,k2,k3) point, as a 3D matrix
for ik=-Nhalf:Nhalf
  k1               = kunit*ik;
  iksft            = ik + Nhalf+1;
  k1_3D(iksft,:,:) = k1;
end
%% k2 component on each (k1,k2,k3) point, as a 3D matrix
for jk=-Nhalf:Nhalf
  k2               = kunit*jk;
  jksft            = jk + Nhalf+1;
  k2_3D(:,jksft,:) = k2;
end
%% k3 component on each (k1,k2,k3) point, as a 3D matrix
for kk=-Nhalf:Nhalf
  k3               = kunit*kk;
  kksft            = kk + Nhalf+1;
  k3_3D(:,:,kksft) = k3;
end

ksq     = k1_3D.^2 +k2_3D.^2 +k3_3D.^2; %% 3D matrix of k^2.

%% 3D k-space fluctuations at z=1000, before applying random seed
Deltacval = interp1(kktab, Dc_zi,  sqrt(ksq), 'spline', 'extrap'); %% 3D matrix 
Deltabval = interp1(kktab, Db_zi,  sqrt(ksq), 'spline', 'extrap'); %% 3D matrix
Deltarval = interp1(kktab, Dr_zi,  sqrt(ksq), 'spline', 'extrap'); %% 3D matrix

Thetacval = interp1(kktab, THc_zi, sqrt(ksq), 'spline', 'extrap'); %% 3D matrix
Thetabval = interp1(kktab, THb_zi, sqrt(ksq), 'spline', 'extrap'); %% 3D matrix
Vcbval    = interp1(kktab, Vcb_zi, sqrt(ksq), 'spline', 'extrap'); %% 3D matrix

DeltaTval = interp1(kktab, DT_zi,  sqrt(ksq), 'spline', 'extrap'); %% 3D matrix

%% Get modes in 3D k-space matrices:
Deltacomval  = interp1(kktab, Deltacom_k, sqrt(ksq), 'spline', 'extrap');
Deltastrval  = interp1(kktab, Deltastr_k, sqrt(ksq), 'spline', 'extrap');
Deltagroval  = interp1(kktab, Deltagro_k, sqrt(ksq), 'spline', 'extrap');
Deltadecval  = interp1(kktab, Deltadec_k, sqrt(ksq), 'spline', 'extrap');

%% Now apply random seed.
%% /sqrt(2) is for distributing P(k) to both real and imaginary
%% See e.g. astro-ph/0506540, equation (60) & (61).
Delta_c_k    = Deltacval          .* gauss3D    /sqrt(2);
Theta_c_k    = Thetacval          .* gauss3D    /sqrt(2);
V_c_k_1      = -i*ai*k1_3D./ksq   .* Theta_c_k          ; %% sqrt(2) included above
V_c_k_2      = -i*ai*k2_3D./ksq   .* Theta_c_k          ; %% sqrt(2) included above
V_c_k_3      = -i*ai*k3_3D./ksq   .* Theta_c_k          ; %% sqrt(2) included above
Delta_b_k    = Deltabval          .* gauss3D    /sqrt(2);
Theta_b_k    = Thetabval          .* gauss3D    /sqrt(2);
V_b_k_1      = -i*ai*k1_3D./ksq   .* Theta_b_k          ; %% sqrt(2) included above
V_b_k_2      = -i*ai*k2_3D./ksq   .* Theta_b_k          ; %% sqrt(2) included above
V_b_k_3      = -i*ai*k3_3D./ksq   .* Theta_b_k          ; %% sqrt(2) included above
Delta_T_k    = DeltaTval          .* gauss3D    /sqrt(2);
Delta_com_k  = Deltacomval        .* gauss3D    /sqrt(2);
Delta_str_k  = Deltastrval        .* gauss3D    /sqrt(2);
Delta_gro_k  = Deltagroval        .* gauss3D    /sqrt(2);
Delta_dec_k  = Deltadecval        .* gauss3D    /sqrt(2);
%% V_cb_k_* = V_c_k_* - V_b_k_*, or we can do the following.
%V_cb_k_1      = -i*k1_3D./sqrt(ksq) .*V_cb_val  .* gauss3D    /sqrt(2);
%V_cb_k_2      = -i*k2_3D./sqrt(ksq) .*V_cb_val  .* gauss3D    /sqrt(2);
%V_cb_k_3      = -i*k3_3D./sqrt(ksq) .*V_cb_val  .* gauss3D    /sqrt(2);
%%----- Accurate initialization, using time derivative of CAMB data --- begin

%% Apply reality condition on face (kk=0), across x axis.
%% for CDM
Delta_c_k  (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_c_k  (1:Nmode, 1:Nc-1, Nc));
Theta_c_k  (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Theta_c_k  (1:Nmode, 1:Nc-1, Nc));
V_c_k_1    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_c_k_1    (1:Nmode, 1:Nc-1, Nc));
V_c_k_2    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_c_k_2    (1:Nmode, 1:Nc-1, Nc));
V_c_k_3    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_c_k_3    (1:Nmode, 1:Nc-1, Nc));
%% for baryon						         
Delta_b_k  (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_b_k  (1:Nmode, 1:Nc-1, Nc));
Theta_b_k  (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Theta_b_k  (1:Nmode, 1:Nc-1, Nc));
V_b_k_1    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_b_k_1    (1:Nmode, 1:Nc-1, Nc));
V_b_k_2    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_b_k_2    (1:Nmode, 1:Nc-1, Nc));
V_b_k_3    (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(V_b_k_3    (1:Nmode, 1:Nc-1, Nc));
Delta_T_k  (Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_T_k  (1:Nmode, 1:Nc-1, Nc));
%% for modes
Delta_com_k(Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_com_k(1:Nmode, 1:Nc-1, Nc));
Delta_str_k(Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_str_k(1:Nmode, 1:Nc-1, Nc));
Delta_gro_k(Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_gro_k(1:Nmode, 1:Nc-1, Nc));
Delta_dec_k(Nmode:-1:1, Nmode:-1:Nc+1, Nc) = conj(Delta_dec_k(1:Nmode, 1:Nc-1, Nc));

%% Apply reality condition on axis (jj=0, kk=0), along x axis.
%% for CDM
Delta_c_k  (Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_c_k  (1:Nc-1, Nc, Nc));
Theta_c_k  (Nmode:-1:Nc+1, Nc, Nc) = conj(Theta_c_k  (1:Nc-1, Nc, Nc));
V_c_k_1    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_c_k_1    (1:Nc-1, Nc, Nc));
V_c_k_2    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_c_k_2    (1:Nc-1, Nc, Nc));
V_c_k_3    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_c_k_3    (1:Nc-1, Nc, Nc));
%% for baryon					      	 
Delta_b_k  (Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_b_k  (1:Nc-1, Nc, Nc));
Theta_b_k  (Nmode:-1:Nc+1, Nc, Nc) = conj(Theta_b_k  (1:Nc-1, Nc, Nc));
V_b_k_1    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_b_k_1    (1:Nc-1, Nc, Nc));
V_b_k_2    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_b_k_2    (1:Nc-1, Nc, Nc));
V_b_k_3    (Nmode:-1:Nc+1, Nc, Nc) = conj(V_b_k_3    (1:Nc-1, Nc, Nc));
Delta_T_k  (Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_T_k  (1:Nc-1, Nc, Nc));
%% for Delta_{-} coefficients			      	 
Delta_com_k(Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_com_k(1:Nc-1, Nc, Nc));
Delta_str_k(Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_str_k(1:Nc-1, Nc, Nc));
Delta_gro_k(Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_gro_k(1:Nc-1, Nc, Nc));
Delta_dec_k(Nmode:-1:Nc+1, Nc, Nc) = conj(Delta_dec_k(1:Nc-1, Nc, Nc));

%% Now for the lower half apply reality condition
%% for CDM
Delta_c_k  (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_c_k  (1:Nmode, 1:Nmode, 1:Nc-1));
Theta_c_k  (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Theta_c_k  (1:Nmode, 1:Nmode, 1:Nc-1));
V_c_k_1    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_c_k_1    (1:Nmode, 1:Nmode, 1:Nc-1));
V_c_k_2    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_c_k_2    (1:Nmode, 1:Nmode, 1:Nc-1));
V_c_k_3    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_c_k_3    (1:Nmode, 1:Nmode, 1:Nc-1));
%% for baryon
Delta_b_k  (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_b_k  (1:Nmode, 1:Nmode, 1:Nc-1));
Theta_b_k  (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Theta_b_k  (1:Nmode, 1:Nmode, 1:Nc-1));
V_b_k_1    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_b_k_1    (1:Nmode, 1:Nmode, 1:Nc-1));
V_b_k_2    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_b_k_2    (1:Nmode, 1:Nmode, 1:Nc-1));
V_b_k_3    (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(V_b_k_3    (1:Nmode, 1:Nmode, 1:Nc-1));
Delta_T_k  (Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_T_k  (1:Nmode, 1:Nmode, 1:Nc-1));
%% for Delta_{-} coefficients
Delta_com_k(Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_com_k(1:Nmode, 1:Nmode, 1:Nc-1));
Delta_str_k(Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_str_k(1:Nmode, 1:Nmode, 1:Nc-1));
Delta_gro_k(Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_gro_k(1:Nmode, 1:Nmode, 1:Nc-1));
Delta_dec_k(Nmode:-1:1, Nmode:-1:1, Nmode:-1:Nc+1) = conj(Delta_dec_k(1:Nmode, 1:Nmode, 1:Nc-1));

%% Nullify the monopole term
Delta_c_k  (Nc,Nc,Nc) = complex(0);
Theta_c_k  (Nc,Nc,Nc) = complex(0);
V_c_k_1    (Nc,Nc,Nc) = complex(0);
V_c_k_2    (Nc,Nc,Nc) = complex(0);
V_c_k_3    (Nc,Nc,Nc) = complex(0);
Delta_b_k  (Nc,Nc,Nc) = complex(0);
Theta_b_k  (Nc,Nc,Nc) = complex(0);
V_b_k_1    (Nc,Nc,Nc) = complex(0);
V_b_k_2    (Nc,Nc,Nc) = complex(0);
V_b_k_3    (Nc,Nc,Nc) = complex(0);
Delta_T_k  (Nc,Nc,Nc) = complex(0);
Delta_com_k(Nc,Nc,Nc) = complex(0);
Delta_str_k(Nc,Nc,Nc) = complex(0);
Delta_gro_k(Nc,Nc,Nc) = complex(0);
Delta_dec_k(Nc,Nc,Nc) = complex(0);


%% multiply forgotten coefficient.
%% sqrt(P(k)) has unit of sqrt(volume), so to have dimensionless Delta_c_k
%% divide below by sqrt(Vbox). (See also Coles&Luccin for how Delta_k is defined)
%% 3D fft of matlab defined as Ak(k_n) = Sigma_j A(x_j) exp(-i * k_n dot x_j), and
%% 3D ifft of matlab defined as A(x_j) = 1/N^3 Sigma_n Ak(k_n) exp(i * k_n dot x_j)
%% while typical cosmology(e.g. Coles&Luccin) uses
%% Ak(k_n) = 1/V   int d^3x A(x_j) exp(-i * k_n dot x_j)
%%         = 1/N^3 Sigma_j  A(x_j) exp(-i * k_n dot x_j)
%% So, Ak(k_n, matlab) * N^3 = Ak(k_n, cosmology)
%% See e.g. astro-ph/0506540, equation (60) & (61). (1/sqrt(2) included above)
Delta_c_k   = Delta_c_k   /sqrt(Vbox) *Nmode^3;  %% unitless
Theta_c_k   = Theta_c_k   /sqrt(Vbox) *Nmode^3;  %% 1/Myr   
V_c_k_1     = V_c_k_1     /sqrt(Vbox) *Nmode^3;  %% Mpc/Myr 
V_c_k_2     = V_c_k_2     /sqrt(Vbox) *Nmode^3;  
V_c_k_3     = V_c_k_3     /sqrt(Vbox) *Nmode^3;
Delta_b_k   = Delta_b_k   /sqrt(Vbox) *Nmode^3;  %% unitless
Theta_b_k   = Theta_b_k   /sqrt(Vbox) *Nmode^3;  %% 1/Myr   
V_b_k_1     = V_b_k_1     /sqrt(Vbox) *Nmode^3;  %% Mpc/Myr 
V_b_k_2     = V_b_k_2     /sqrt(Vbox) *Nmode^3;
V_b_k_3     = V_b_k_3     /sqrt(Vbox) *Nmode^3;
Delta_T_k   = Delta_T_k   /sqrt(Vbox) *Nmode^3;  %% unitless
Delta_com_k = Delta_com_k /sqrt(Vbox) *Nmode^3;
Delta_str_k = Delta_str_k /sqrt(Vbox) *Nmode^3;
Delta_gro_k = Delta_gro_k /sqrt(Vbox) *Nmode^3;
Delta_dec_k = Delta_dec_k /sqrt(Vbox) *Nmode^3;


%% When taking ifft, first shift k-space matrix to "default" matlab k-space one.
%% From the beginning I constructed monopole-centered k-space matrix.
Delta_c_k   = ifftshift(Delta_c_k    );
Theta_c_k   = ifftshift(Theta_c_k    );
V_c_k_1     = ifftshift(V_c_k_1      );
V_c_k_2     = ifftshift(V_c_k_2      );
V_c_k_3     = ifftshift(V_c_k_3      );
Delta_b_k   = ifftshift(Delta_b_k    );
Theta_b_k   = ifftshift(Theta_b_k    );
V_b_k_1     = ifftshift(V_b_k_1      );
V_b_k_2     = ifftshift(V_b_k_2      );
V_b_k_3     = ifftshift(V_b_k_3      );
Delta_T_k   = ifftshift(Delta_T_k    );
Delta_com_k = ifftshift(Delta_com_k );
Delta_str_k = ifftshift(Delta_str_k );
Delta_gro_k = ifftshift(Delta_gro_k  );
Delta_dec_k = ifftshift(Delta_dec_k);

%% Now, the initial real-space fields.
Delta_c       = real(ifftn(Delta_c_k    ));
Theta_c       = real(ifftn(Theta_c_k    ));
V_c_1         = real(ifftn(V_c_k_1      ));
V_c_2         = real(ifftn(V_c_k_2      ));
V_c_3         = real(ifftn(V_c_k_3      ));
Delta_b       = real(ifftn(Delta_b_k    ));
Theta_b       = real(ifftn(Theta_b_k    ));
V_b_1         = real(ifftn(V_b_k_1      ));
V_b_2         = real(ifftn(V_b_k_2      ));
V_b_3         = real(ifftn(V_b_k_3      ));
Delta_T       = real(ifftn(Delta_T_k    ));
Deltacom      = real(ifftn(Delta_com_k ));
Deltastr      =	real(ifftn(Delta_str_k ));
Deltagro      =	real(ifftn(Delta_gro_k  ));
Deltadec      =	real(ifftn(Delta_dec_k));
V_cb_1        = V_c_1 - V_b_1;
V_cb_2        = V_c_2 - V_b_2;
V_cb_3        = V_c_3 - V_b_3;
Vcb           = sqrt(V_cb_1.^2+V_cb_2.^2+V_cb_3.^2);

%% BCCOMICS needs large-scale monopole values at output redshift -- begin
%% First at z=1000
%% See last comment in Integrate_evolODE.m why we save in two formats.
if matlabflag
  %% to matlab binary format
  save([outputdir '/Dc3D.matbin'],   'Delta_c', '-v6'); 
  save([outputdir '/Db3D.matbin'],   'Delta_b', '-v6'); 
  save([outputdir '/THc3D.matbin'],  'Theta_c', '-v6'); 
  save([outputdir '/THb3D.matbin'],  'Theta_b', '-v6'); 
  save([outputdir '/DT.matbin'],     'Delta_T', '-v6');

  save([outputdir '/V_cb_1.matbin'], 'V_cb_1',  '-v6'); 
  save([outputdir '/V_cb_2.matbin'], 'V_cb_2',  '-v6'); 
  save([outputdir '/V_cb_3.matbin'], 'V_cb_3',  '-v6'); 
  save([outputdir '/V_c_1.matbin'],  'V_c_1',   '-v6'); 
  save([outputdir '/V_c_2.matbin'],  'V_c_2',   '-v6'); 
  save([outputdir '/V_c_3.matbin'],  'V_c_3',   '-v6'); 

  %% to hdf5 format
  save([outputdir '/Dc3D.h5'],   'Delta_c', '-v7.3'); 
  save([outputdir '/Db3D.h5'],   'Delta_b', '-v7.3'); 
  save([outputdir '/THc3D.h5'],  'Theta_c', '-v7.3'); 
  save([outputdir '/THb3D.h5'],  'Theta_b', '-v7.3'); 
  save([outputdir '/DT.h5'],     'Delta_T', '-v7.3');

  save([outputdir '/V_cb_1.h5'], 'V_cb_1',  '-v7.3'); 
  save([outputdir '/V_cb_2.h5'], 'V_cb_2',  '-v7.3'); 
  save([outputdir '/V_cb_3.h5'], 'V_cb_3',  '-v7.3'); 
  save([outputdir '/V_c_1.h5'],  'V_c_1',   '-v7.3'); 
  save([outputdir '/V_c_2.h5'],  'V_c_2',   '-v7.3'); 
  save([outputdir '/V_c_3.h5'],  'V_c_3',   '-v7.3'); 
else
  %% to matlab binary format
  save('-mat-binary', [outputdir '/Dc3D.matbin'],   'Delta_c'); 
  save('-mat-binary', [outputdir '/Db3D.matbin'],   'Delta_b'); 
  save('-mat-binary', [outputdir '/THc3D.matbin'],  'Theta_c'); 
  save('-mat-binary', [outputdir '/THb3D.matbin'],  'Theta_b'); 
  save('-mat-binary', [outputdir '/DT.matbin'],     'Delta_T');

  save('-mat-binary', [outputdir '/V_cb_1.matbin'], 'V_cb_1'); 
  save('-mat-binary', [outputdir '/V_cb_2.matbin'], 'V_cb_2'); 
  save('-mat-binary', [outputdir '/V_cb_3.matbin'], 'V_cb_3'); 
  save('-mat-binary', [outputdir '/V_c_1.matbin'],  'V_c_1' ); 
  save('-mat-binary', [outputdir '/V_c_2.matbin'],  'V_c_2' ); 
  save('-mat-binary', [outputdir '/V_c_3.matbin'],  'V_c_3' ); 

  %% to hdf5 format
  save('-hdf5', [outputdir '/Dc3D.h5'],   'Delta_c'); 
  save('-hdf5', [outputdir '/Db3D.h5'],   'Delta_b'); 
  save('-hdf5', [outputdir '/THc3D.h5'],  'Theta_c'); 
  save('-hdf5', [outputdir '/THb3D.h5'],  'Theta_b'); 
  save('-hdf5', [outputdir '/DT.h5'],     'Delta_T');

  save('-hdf5', [outputdir '/V_cb_1.h5'], 'V_cb_1'); 
  save('-hdf5', [outputdir '/V_cb_2.h5'], 'V_cb_2'); 
  save('-hdf5', [outputdir '/V_cb_3.h5'], 'V_cb_3'); 
  save('-hdf5', [outputdir '/V_c_1.h5'],  'V_c_1' ); 
  save('-hdf5', [outputdir '/V_c_2.h5'],  'V_c_2' ); 
  save('-hdf5', [outputdir '/V_c_3.h5'],  'V_c_3' ); 
end

if plotflag
  Plot_slices;  %%==== script ==================
end
