%% Generate 3D real-space maps at z=zzend

V_cb_1_azend = V_cb_1 * (azbegin/azend);
V_cb_2_azend = V_cb_2 * (azbegin/azend);
V_cb_3_azend = V_cb_3 * (azbegin/azend);

a           = azend;
Get_D_dDda;  %%==== script: should always be preceeded by scale factor a.
Dc3D_azend  = Deltagro*Dpg + Deltadec*Dpd +fb*(Deltacom + Deltastr*Dms);
Db3D_azend  = Deltagro*Dpg + Deltadec*Dpd -fc*(Deltacom + Deltastr*Dms);
THc3D_azend = -aH*(Deltagro*dDpg_da + Deltadec*dDpd_da) -aH*fb*Deltastr*dDms_da;
THb3D_azend = -aH*(Deltagro*dDpg_da + Deltadec*dDpd_da) +aH*fc*Deltastr*dDms_da;

%% For DT, first get DTA(at a=0.01) and DTB(at a=0.1), and use
%% the empirical fit for 4Mpc cell to obtain DeltaT at azend. 
%% This also prepares for getting DeltaT at any a.
Get_DeltaT_fit; %%==== script ==================

if matlabflag
  save([outputdir '/V_cb_1_azend.dat'], 'V_cb_1_azend', '-v6'); % @ azend
  save([outputdir '/V_cb_2_azend.dat'], 'V_cb_2_azend', '-v6'); % @ azend
  save([outputdir '/V_cb_3_azend.dat'], 'V_cb_3_azend', '-v6'); % @ azend
  save([outputdir '/Dc3D_azend.dat'],   'Dc3D_azend',   '-v6'); % @ azend
  save([outputdir '/Db3D_azend.dat'],   'Db3D_azend',   '-v6'); % @ azend
  save([outputdir '/THc3D_azend.dat'],  'THc3D_azend',  '-v6'); % @ azend
  save([outputdir '/THb3D_azend.dat'],  'THb3D_azend',  '-v6'); % @ azend
  save([outputdir '/DT_azend.dat'],     'DT3D_azend',   '-v6'); % @ azend
else
  save('-mat-binary', [outputdir '/V_cb_1_azend.dat'], 'V_cb_1_azend'); % @ azend
  save('-mat-binary', [outputdir '/V_cb_2_azend.dat'], 'V_cb_2_azend'); % @ azend
  save('-mat-binary', [outputdir '/V_cb_3_azend.dat'], 'V_cb_3_azend'); % @ azend
  save('-mat-binary', [outputdir '/Dc3D_azend.dat'],   'Dc3D_azend'  ); % @ azend
  save('-mat-binary', [outputdir '/Db3D_azend.dat'],   'Db3D_azend'  ); % @ azend
  save('-mat-binary', [outputdir '/THc3D_azend.dat'],  'THc3D_azend' ); % @ azend
  save('-mat-binary', [outputdir '/THb3D_azend.dat'],  'THb3D_azend' ); % @ azend
  save('-mat-binary', [outputdir '/DT_azend.dat'],     'DT3D_azend'  ); % @ azend
end

if plotflag
  ifig = ifig+1;
  figure(ifig);
  if ~matflag
    disp('Be patient, hist2 plot takes a few minutes');
  end
  hist2d([Dc3D_azend(:),THc3D_azend(:)],100,100);
  daspect([1 1 1]);
  view(2);
end
