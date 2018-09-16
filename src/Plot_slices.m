%% Plots some quantities on a 2D slice

Zc   = reshape(Delta_c(1,:,:),Nmode,Nmode);
Zthc = reshape(Theta_c(1,:,:),Nmode,Nmode);
Zb   = reshape(Delta_b(1,:,:),Nmode,Nmode);
ZT   = reshape(Delta_T(1,:,:),Nmode,Nmode);

ifig = ifig+1;
figure(ifig);
imagesc(flipud(rot90(Zc  ))), colorbar;
title('\Delta_{c}');

%% Should look almost opposite of Delta_c plot, because Theta_c is almost
%% proportional to -Delta_c. See Fig.2(d) of A16.
ifig = ifig+1;
figure(ifig);
imagesc(flipud(rot90(Zthc))), colorbar;
title('\Theta_{c} (Myr^{-1})');

ifig = ifig+1;
figure(ifig);
imagesc(flipud(rot90(Zb  ))), colorbar;
title('\Delta_{b}');

ifig = ifig+1;
figure(ifig);
imagesc(flipud(rot90(ZT  ))), colorbar;
title('\Delta_{T}');
