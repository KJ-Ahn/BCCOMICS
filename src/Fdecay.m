%% Eq. 27 of A16, after calculating dlnH/da and Omega_matter

function dFda = Fdecay(a,F)
  global Om0 Omr0 OmLambda0;

  dFda    = zeros(2,1);
  dFda(1) = F(2);
  dFda(2) = (3*Om0*a + 4*Omr0)/(2*a*(Om0*a+Omr0+OmLambda0*a^4))*F(2) + (-9*Omr0+3*OmLambda0*a^4)/(4*a^2*(Om0*a+Omr0+OmLambda0*a^4))*F(1);
end

