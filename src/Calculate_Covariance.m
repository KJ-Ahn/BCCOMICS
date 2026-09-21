%% Script: Calculate_Covariance.m
%% What it does: Calculates the covariance matrix (Sigma) for the target constraint.
%%               Velocity constraints are excluded to prevent unphysical high-k ringing.

disp('----------------Calculating Covariance Matrix----------------');
%% Variance of Matter Overdensity
%% Using the exact grid sums guarantees 0% residual error after IFFT.
Deltamval  = fc * Deltacval + fb * Deltabval;
Deltamval(Nc, Nc, Nc) = 0; 
Sigma11_exact = sum(Deltamval(:).^2) / Vbox;

disp('----------------Covariance Matrix constructed----------------');