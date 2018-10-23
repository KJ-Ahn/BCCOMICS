%% integrant d(a_local)/d(t*H_i), where H_i is NOT the local initial Hubble but
%% is the global initial Hubble.
function daloc_dtHi = fdadt(tHi, aloc)
  global H_i H_l_i Om_l_i Omr_l_i OmLambda_l_i OmK_l_i aloci;
  Hratio = H_l_i/H_i;
  
  daloc_dtHi = Hratio * aloci* sqrt(Om_l_i/(aloc/aloci) + Omr_l_i/(aloc/aloci)^2 + OmLambda_l_i*(aloc/aloci)^2 + OmK_l_i);
end