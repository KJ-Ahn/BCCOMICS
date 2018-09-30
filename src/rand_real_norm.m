%% Function that obtains random phase and then applying reality condition
%% to realize k-space deltas. Ak and Bk are (N,N,N) complex matrices.
%% ramp, rphs are 1D arrays of N*N*Nc elements. Nc=N/2+1 with even N.
%% For reality condition, the along-axis-largest k modes 
%% (e.g. kx=Nc*k_unit) do not have matching pairs if 
%% N is even. Just leave them out.
%% Also multipy left-out normalization
function Bk = rand_real_norm(Ak,N,Nc,ramp,rphs,Vbox)
  Bk           = zeros(N,N,N);


  %%%%%% All points except three leftmost faces  ------------------ begin
  %% Apply (identical) random seed
  Ak(:,:,1:Nc) = Ak(:,:,1:Nc).*reshape(ramp(:).*exp(i*rphs(:)),N,N,Nc);
  Bk           = Ak;
  %% Apply reality condition on face (kk=0), across x axis.
  Bk (N:-1:2, N:-1:Nc+1, Nc)     = conj(Ak (2:N, 2:Nc-1, Nc));
  %% Apply reality condition on axis (jj=0, kk=0), along x axis.
  Bk (N:-1:Nc+1, Nc, Nc)         = conj(Ak (2:Nc-1, Nc, Nc));

  %% Now for the lower half apply reality condition but excluding kk=0 plane
  Bk (N:-1:2, N:-1:2, N:-1:Nc+1) = conj(Ak (2:N, 2:N, 2:Nc-1));
  %% Nullify monopole terms. 
  Bk (Nc,Nc,Nc) = complex(0);
  %%%%%% All points except three leftmost faces  ------------------ end


  %% --Apply another reality condition for even # sampled k values.
  %% --Because ik=-Nhalf_H, etc. is the symmetric point, it is guaranteed to
  %% --have the conjugate relation as for ii=0, jj=0, kk=0 face-on cases.

  %%%%%% Three left-out faces  -------------------------- begin
  %%%% On face  (ii=-N/2), without leftmost edges
  %% Apply reality condition across projected y axis.
  Bk (1, N:-1:2, N:-1:Nc+1)     = conj(Ak (1, 2:N, 2:Nc-1));
  %% Apply reality condition on axis (kk=0), along projected y axis.
  Bk (1, N:-1:Nc+1, Nc)         = conj(Ak (1, 2:Nc-1, Nc));

  %%%% On face  (jj=-N/2), without leftmost edges
  %% Apply reality condition across projected y axis.
  Bk (N:-1:2, 1, N:-1:Nc+1)     = conj(Ak (2:N, 1, 2:Nc-1));
  %% Apply reality condition on axis (kk=0), along projected y axis.
  Bk (N:-1:Nc+1, 1, Nc)         = conj(Ak (2:Nc-1, 1, Nc));

  %%%% On face  (kk=-N/2), without leftmost edges
  %% Apply reality condition across projected x axis.
  Bk (N:-1:2, N:-1:Nc+1, 1)     = conj(Ak (2:N, 2:Nc-1, 1));
  %% Apply reality condition on axis (jj=0), along projected x axis.
  Bk (N:-1:Nc+1, Nc, 1)         = conj(Ak (2:Nc-1, Nc, 1));

  %%%% On leftmost edges
  %% Apply reality condition on axis (jj=-N/2), along x-edge.
  Bk (N:-1:Nc+1, 1, 1)         = conj(Ak (2:Nc-1, 1, 1));

  %% Apply reality condition on axis (ii=-N/2), along y-edge.
  Bk (1, N:-1:Nc+1, 1)         = conj(Ak (1, 2:Nc-1, 1));

  %% Apply reality condition on axis (ii=-N/2), along z-edge.
  Bk (1, 1, N:-1:Nc+1)         = conj(Ak (1, 1, 2:Nc-1));

  %%%% On symmetry point (it should be a real number)
  Bk(1,1,1)                    = abs(Bk(1,1,1));
  %%%%%% Three left-out faces  -------------------------- end

  %%%%%% Normalization
  Bk = Bk /sqrt(Vbox)*N^3;
  
end
