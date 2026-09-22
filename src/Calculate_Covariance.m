%% Script: Calculate_Covariance.m
%% What it does: Calculates the covariance matrix (Sigma) for the target density constraint.
%%
%% [METHODOLOGICAL NOTE FOR PEER REVIEW]
%% Why is the velocity (V_cb) constraint excluded from the Hoffman-Ribak (HR) filter?
%% If V_cb is constrained via the HR filter, its inherent cross-correlation with the 
%% streaming mode and baryon velocity will alter the local density field. This breaks 
%% the strict 'identical density background' condition required for our controlled 
%% experiment. Therefore, density is constrained via HR, while the target V_cb is 
%% achieved via a uniform box-wide offset (k=0 component) in Apply_CR_filter.m, 
%% ensuring the underlying scalar fields remain bit-for-bit identical across V_cb cases.

disp('----------------Calculating Covariance Matrix----------------');
%% Variance of Matter Overdensity
%% Using the exact grid sums guarantees 0% residual error after IFFT.
Deltamval  = fc * Deltacval + fb * Deltabval;
Deltamval(Nc, Nc, Nc) = 0; 
Sigma11_exact = sum(Deltamval(:).^2) / Vbox;

disp('----------------Covariance Matrix constructed----------------');
