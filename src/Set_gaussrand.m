%% Create/Use/Shuffle/Add random numbers for initial condition
%% 'subgaussseed512.matbin' is used as the seed for any Ncell_p<=512 cases,
%% when oldseedflag is true. For Ncell_p>512 cases, if oldseedflag is still true,
%% 'subgaussseed512.matbin' is used as part of the seed and missing seeds are
%% generated attached.
%% Seeds are arranged such that resolution studies can be done.

file512seed = [setupdir '/subgaussseed512.matbin'];
fileNseed   = [setupdir '/subgaussseed' num2str(Nmode_p) '.matbin'];
if oldseedflag  %% use preexisting seed
  if ~exist(file512seed)
    disp(['You let oldseedflag true, but ' setupdir '/subgaussseed512.matbin does not exist.']);
    returnflag=true;
    return;
  end
  if (Nmode_p <= 512)  %% use preexisting seed for Nmode_p <= 512
    if matlabflag
      load(file512seed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', file512seed, 'randamp', 'randphs');
    end
    clear file512seed;  %% save memory 
    istart = (512-Nmode_p)/2 + 1;  %% x & y starting index of 512*512*257 seed to use
    iend   = istart + Nmode_p -1;  %% x & y ending index of 512*512*257 seed to use
    randamp_ = randamp(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
    randphs_ = randphs(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
    clear randamp randphs;  %% save memory 
    randamp = reshape(randamp_, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
    randphs = reshape(randphs_, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
    clear randamp_ randphs_;  %% save memory
  elseif exist(fileNseed)  %% use preexisting seed, and also generate & attach higher-k seeds.
    if matlabflag
      load(fileNseed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', fileNseed, 'randamp', 'randphs');
    end
    clear fileNseed;
  elseif ~exist(fileNseed)  %% use preexisting seed, and also generate & attach higher-k seeds.
    if matlabflag
      load(file512seed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', file512seed, 'randamp', 'randphs');
    end
    clear file512seed;
    Nmissing = Nmode_p*Nmode_p*Nc_p - 512*512*257; %% # of missing modes
    randamp_  = zeros(Nmode_p, Nmode_p, Nc_p);   %% seed to fill in
    randphs_  = zeros(Nmode_p, Nmode_p, Nc_p);   %% seed to fill in
    randamp__ = raylrnd(1,      [Nmissing,1]);   %% new seed (amplitude)
    randphs__ = unifrnd(0.2*pi, [Nmissing,1]);   %% new seed (phase)
    istart = (Nmode_p-512)/2 + 1;   %% x & y starting index to host 512*512*257 seed
    iend   = istart    + 512 - 1;   %% x & y ending index to host 512*512*257 seed
    izstart = Nc_p-257+1;           %% z starting index to host 512*512*257 seed

    %% host subgaussseed512.matbin seed
    randamp_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randamp, 512, 512, 257);
    randphs_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randphs, 512, 512, 257);
    
    %% attach - 1 (bottom)
    is    = 1;                           %% starting index in 1D seed array
    Ntemp = Nmode_p*Nmode_p*(Nc_p-257);  %% # of modes being filled in
    ie    = is + Ntemp-1;                %% ending index in 1D seed array  
    randamp_(1:Nmode_p, 1:Nmode_p, 1:Nc_p-257)        = reshape(randamp__(is:ie), Nmode_p, Nmode_p, Nc_p-257);  %% fill in part of randamp_
    randphs_(1:Nmode_p, 1:Nmode_p, 1:Nc_p-257)        = reshape(randphs__(is:ie), Nmode_p, Nmode_p, Nc_p-257);  %% fill in part of randphs_

    %% attach - 2 (top left)
    is    = ie + 1;                  %% starting index in 1D seed array
    Ntemp = (istart-1)*Nmode_p*257;  %% # of modes being filled in     
    ie    = is + Ntemp-1;            %% ending index in 1D seed array  
    randamp_(1:istart-1, 1:Nmode_p, izstart:Nc_p)     = reshape(randamp__(is:ie), istart-1, Nmode_p, 257);  %% fill in part of randamp_
    randphs_(1:istart-1, 1:Nmode_p, izstart:Nc_p)     = reshape(randphs__(is:ie), istart-1, Nmode_p, 257);  %% fill in part of randphs_

    %% attach - 3 (top right)
    is    = ie + 1;                     %% starting index in 1D seed array
    Ntemp = (Nmode_p-iend)*Nmode_p*257; %% # of modes being filled in     
    ie    = is + Ntemp-1;               %% ending index in 1D seed array  
    randamp_(iend+1:Nmode_p, 1:Nmode_p, izstart:Nc_p) = reshape(randamp__(is:ie), Nmode_p-iend, Nmode_p, 257);  %% fill in part of randamp_
    randphs_(iend+1:Nmode_p, 1:Nmode_p, izstart:Nc_p) = reshape(randphs__(is:ie), Nmode_p-iend, Nmode_p, 257);  %% fill in part of randphs_

    %% attach - 4 (top middle front)
    is    = ie + 1;                         %% starting index in 1D seed array
    Ntemp = (iend-istart+1)*(istart-1)*257; %% # of modes being filled in     
    ie    = is + Ntemp-1;                   %% ending index in 1D seed array  
    randamp_(istart:iend, 1:istart-1, izstart:Nc_p)   = reshape(randamp__(is:ie), iend-istart+1, istart-1, 257);  %% fill in part of randamp_
    randphs_(istart:iend, 1:istart-1, izstart:Nc_p)   = reshape(randphs__(is:ie), iend-istart+1, istart-1, 257);  %% fill in part of randphs_

    %% attach - 5 (top middle back)
    is    = ie + 1;                             %% starting index in 1D seed array
    Ntemp = (iend-istart+1)*(Nmode_p-iend)*257; %% # of modes being filled in     
    ie    = is + Ntemp-1;                       %% ending index in 1D seed array  
    randamp_(istart:iend, iend+1:Nmode_p, izstart:Nc_p)   = reshape(randamp__(is:ie), iend-istart+1, Nmode_p-iend, 257);  %% fill in part of randamp_
    randphs_(istart:iend, iend+1:Nmode_p, izstart:Nc_p)   = reshape(randphs__(is:ie), iend-istart+1, Nmode_p-iend, 257);  %% fill in part of randphs_

    clear randamp randphs randamp__ randphs__;
    randamp = reshape(randamp_, Nmode_p*Nmode_p*Nc_p, 1);
    randphs = reshape(randphs_, Nmode_p*Nmode_p*Nc_p, 1);
    clear randamp_ randphs_;
    
  end
else  %% generate new seed
  randamp = raylrnd(1,      [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (amplitude)
  randphs = unifrnd(0.2*pi, [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (phase)
end

