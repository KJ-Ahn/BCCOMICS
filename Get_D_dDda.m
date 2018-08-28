%%%% Script: Calculates (scale factor) ~ (growth factor) tables.
%%           a is the input.
%%           Needs precalculated tables from Get_growth.m
%% 
%% Author: Kyungjin Ahn
%%
Hz      = H0*sqrt(Om0*a^-3+Omr0*a^-4 + OmLambda0);
aH      = a*Hz;
Dpg     = interp1(azz, Dplus_grow   , a, 'spline');
Dpd     = interp1(azz, Dplus_decay  , a, 'spline');
Dms     = interp1(azz, Dminus_stream, a, 'spline');
dDpg_da = interp1(azz, dDplus_grow_da   , a, 'spline');
dDpd_da = interp1(azz, dDplus_decay_da  , a, 'spline');
dDms_da = interp1(azz, dDminus_stream_da, a, 'spline');
