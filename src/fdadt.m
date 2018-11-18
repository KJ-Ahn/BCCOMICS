%% integrant d(a_local)/d(t*H_i), where H_i is NOT the local initial Hubble but
%% is the global initial Hubble.
function daloc_dtHi = fdadt(tHi, aloc)
  global H_i H_l_i Om_l_i Omr_l_i OmLambda_l_i OmK_l_i aloci;
  Hratio = H_l_i/H_i;
  
  inside_sqrt = Om_l_i/(aloc/aloci) + Omr_l_i/(aloc/aloci)^2 + OmLambda_l_i*(aloc/aloci)^2 + OmK_l_i;
  if (inside_sqrt>=0)
    daloc_dtHi = Hratio * aloci* sqrt(inside_sqrt);
  else
      %% a_local should decay after turnaround in principle, but let a_local frozen
      %% after turnaround. Let enzo not run beyond turnaround.
    daloc_dtHi = 0;
  end
end