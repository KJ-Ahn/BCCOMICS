%% Initialize temperature at z=zi=1000 using approach by NB
%% DT_zi is the product

z1   = 1000.5;
z2   = 999.5;
if (z1<zi | z2>zi)
  disp('Not calculating radiation difference at right redshift');
end

TFr  = TF_zi(:,4);  %% radiation
Pkr  = PS_wo_TFtab .* TFr.^2; %% in Mpc^3 unit

a1   = 1/(1+z1);
a2   = 1/(1+z2);
TF1  = load([TFstr1 num2str(z1) TFstr2]);
TF2  = load([TFstr1 num2str(z2) TFstr2]);
TF1b = TF1(:,3);  %% baryon
TF2b = TF2(:,3);  %% baryon
TF1r = TF1(:,4);  %% radiation
TF2r = TF2(:,4);  %% radiation

Pk1b = PS_wo_TFtab .* TF1b.^2; %% in Mpc^3 unit
Pk2b = PS_wo_TFtab .* TF2b.^2; %% in Mpc^3 unit
Pk1r = PS_wo_TFtab .* TF1r.^2; %% in Mpc^3 unit
Pk2r = PS_wo_TFtab .* TF2r.^2; %% in Mpc^3 unit
dDelta_b_dt = Hzi*ai*(sqrt(Pk2b).*sign(TF2b)-sqrt(Pk1b).*sign(TF1b))/(a2-a1); %% in Mpc^(3/2) Myr^(-1) unit
dDelta_r_dt = Hzi*ai*(sqrt(Pk2r).*sign(TF2r)-sqrt(Pk1r).*sign(TF1r))/(a2-a1); %% in Mpc^(3/2) Myr^(-1) unit

%% Eq. 13 of Ahn13, after NB
DT_zi  = tgamma/xei*ai^4*Tbzi/Tgammai*(2/3*dDelta_b_dt-1/4*dDelta_r_dt) + Dr_zi*(5/4-Tbzi/Tgammai); %% Delta of gas temperature, in units of Mpc^1.5 
PkT    = DT_zi.^2; %% power spectrum in Mpc^3 unit
