%% Script: Apply_CR_filter.m
%% What it does: Applies the Hoffman-Ribak Constrained Realization filter
%%               to the unconstrained k-space fields and updates real-space grids.
%%               Velocity constraints are applied as a global, divergence-free 
%%               shift to prevent unphysical high-k ringing.

disp('-----------------------------------------------------------------------');
disp(' Applying Constrained Realization (CR) Filter in k-space');
disp('-----------------------------------------------------------------------');

%% 1. Extract unconstrained baseline values at target coordinate ---------- begin
unconst_Dm = fc * Delta_c(icc(1), icc(2), icc(3)) + fb * Delta_b(icc(1), icc(2), icc(3));
delta_C1 = target_Dm - unconst_Dm;
%% 1. Extract unconstrained baseline values at target coordinate ---------- end

%% 2. Construct the CR Bending Field (Delta_g_k) -------------------------- begin
dx_shift = (icc(1) - 1) * Lbox / Nmode;
dy_shift = (icc(2) - 1) * Lbox / Nmode;
dz_shift = (icc(3) - 1) * Lbox / Nmode;
phase_shift = exp(-i * (k1_3D * dx_shift + k2_3D * dy_shift + k3_3D * dz_shift));

%% Construct the k-space correction seed (Density only)
Delta_g_k = (1/sqrt(Vbox)) * (Deltamval * (delta_C1 / Sigma11_exact));
Delta_g_k = Delta_g_k .* phase_shift;
Delta_g_k(Nc, Nc, Nc) = complex(0); %% Nullify monopole
%% 2. Construct the CR Bending Field (Delta_g_k) -------------------------- end

%% 3. Apply Correction and Update Real-Space Arrays ----------------------- begin
disp('----- Inverse FFT blending the constraints -----');
norm_factor = 1/sqrt(Vbox) * Nmode^3;
ksq_safe = ksq;
ksq_safe(Nc, Nc, Nc) = 1;

%% Densities & Modes
Delta_c = Delta_c + real(ifftn(ifftshift( Deltacval .* Delta_g_k * norm_factor )));
Theta_c = Theta_c + real(ifftn(ifftshift( Thetacval .* Delta_g_k * norm_factor )));
Delta_b = Delta_b + real(ifftn(ifftshift( Deltabval .* Delta_g_k * norm_factor )));
Theta_b = Theta_b + real(ifftn(ifftshift( Thetabval .* Delta_g_k * norm_factor )));
Delta_T = Delta_T + real(ifftn(ifftshift( DeltaTval .* Delta_g_k * norm_factor )));

Deltacom = Deltacom + real(ifftn(ifftshift( Deltacomval .* Delta_g_k * norm_factor )));
Deltastr = Deltastr + real(ifftn(ifftshift( Deltastrval .* Delta_g_k * norm_factor )));
Deltagro = Deltagro + real(ifftn(ifftshift( Deltagroval .* Delta_g_k * norm_factor )));
Deltadec = Deltadec + real(ifftn(ifftshift( Deltadecval .* Delta_g_k * norm_factor )));

%% Velocities (Update gravitationally induced infall velocities)
V_c_1 = V_c_1 + real(ifftn(ifftshift( (-i*ai*k1_3D./ksq_safe) .* Thetacval .* Delta_g_k * norm_factor )));
V_c_2 = V_c_2 + real(ifftn(ifftshift( (-i*ai*k2_3D./ksq_safe) .* Thetacval .* Delta_g_k * norm_factor )));
V_c_3 = V_c_3 + real(ifftn(ifftshift( (-i*ai*k3_3D./ksq_safe) .* Thetacval .* Delta_g_k * norm_factor )));

V_b_1 = V_b_1 + real(ifftn(ifftshift( (-i*ai*k1_3D./ksq_safe) .* Thetabval .* Delta_g_k * norm_factor )));
V_b_2 = V_b_2 + real(ifftn(ifftshift( (-i*ai*k2_3D./ksq_safe) .* Thetabval .* Delta_g_k * norm_factor )));
V_b_3 = V_b_3 + real(ifftn(ifftshift( (-i*ai*k3_3D./ksq_safe) .* Thetabval .* Delta_g_k * norm_factor )));
%% 3. Apply Correction and Update Real-Space Arrays ----------------------- end

%% 4. Force background bulk flow (V_cb) via Global Divergence-free Shift -- begin
%% Calculate intermediate V_cb
V_cb_1_temp = V_c_1 - V_b_1;
V_cb_2_temp = V_c_2 - V_b_2;
V_cb_3_temp = V_c_3 - V_b_3;

%% Calculate the exact offset required to match the target at the center
offset_1 = target_Vcb_1 - V_cb_1_temp(icc(1), icc(2), icc(3));
offset_2 = target_Vcb_2 - V_cb_2_temp(icc(1), icc(2), icc(3));
offset_3 = target_Vcb_3 - V_cb_3_temp(icc(1), icc(2), icc(3));

%% Apply the offset universally to the CDM velocity field.
%% Because the offset is a constant (k=0 mode), it is perfectly divergence-free
%% and does not violate the continuity equation (Theta remains unchanged).
V_c_1 = V_c_1 + offset_1;
V_c_2 = V_c_2 + offset_2;
V_c_3 = V_c_3 + offset_3;

%% Recalculate Final V_cb
V_cb_1 = V_c_1 - V_b_1;
V_cb_2 = V_c_2 - V_b_2;
V_cb_3 = V_c_3 - V_b_3;
Vcb = sqrt(V_cb_1.^2 + V_cb_2.^2 + V_cb_3.^2);
%% 4. Force background bulk flow (V_cb) via Global Divergence-free Shift -- end

%% 5. Verification Output ------------------------------------------------- begin
check_Dm = fc * Delta_c(icc(1), icc(2), icc(3)) + fb * Delta_b(icc(1), icc(2), icc(3));
check_V1 = V_cb_1(icc(1), icc(2), icc(3)) * MpcMyr_2_kms;
check_V2 = V_cb_2(icc(1), icc(2), icc(3)) * MpcMyr_2_kms;
check_V3 = V_cb_3(icc(1), icc(2), icc(3)) * MpcMyr_2_kms;

disp('-----------------------------------------------------------------------');
disp(' CR Filter Verification (Must match Targets exactly):');
disp([' Final Delta_m: ' num2str(check_Dm)]);
disp([' Final V_cb [x, y, z]: [' num2str(check_V1) '  ' num2str(check_V2) '  ' num2str(check_V3) '] km/s']);
disp('-----------------------------------------------------------------------');
%% 5. Verification Output ------------------------------------------------- end
