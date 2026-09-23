%% Script: Calculate_Covariance.m
%% What it does: Calculates the covariance matrix (Sigma) for the target density constraint.
%%
%% Why is V_cb (= V_c - V_b) not constrained through the Hoffman-Ribak (HR) filter?
%% By isotropy, the zero-lag cross-covariance between a scalar field and a
%% vector field vanishes. An HR constraint on V_cb at the target cell would
%% therefore leave the scalar values of that cell (Delta, Theta, Delta_T, modes)
%% unchanged, but would modify the scalar fields around it. Instead, V_cb is set
%% by a uniform (k=0) offset of the CDM velocity field in Apply_CR_filter.m. This
%% reproduces the HR values at the target cell, and in addition keeps every
%% scalar field of the box bit-for-bit identical across V_cb cases. The price is
%% a non-zero box-averaged V_cb, which is why the statistics of the
%% unconstrained realization are recorded before the CR.

disp('----------------Calculating Covariance Matrix----------------');
%% Variance of Matter Overdensity
%% Using the exact grid sums guarantees 0% residual error after IFFT.
Deltamval  = fc * Deltacval + fb * Deltabval;
Deltamval(Nc, Nc, Nc) = 0; 
Sigma11_exact = sum(Deltamval(:).^2) / Vbox;

disp('----------------Covariance Matrix constructed----------------');
