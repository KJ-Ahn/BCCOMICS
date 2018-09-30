%% Create/Use/Shuffle/Add random numbers for initial condition.
%% First search for fileNseed and use it.
%% In case fileNseed doesnt exist but file512seed does,
%% 'subgaussseed512.matbin' is used as the seed for any Ncell_p<=512 cases,
%% when oldseedflag is true. For Ncell_p>512 cases, if oldseedflag is still true,
%% 'subgaussseed512.matbin' is used as part of the seed and missing seeds are
%% generated and attached.
%% Seeds are arranged such that resolution studies can be done.

if oldseedflag
  record512companionflag = false;
end

file512seed = [setupdir '/subgaussseed512.matbin'];
fileNseed   = [setupdir '/subgaussseed' num2str(Nmode_p) '.matbin'];
if oldseedflag  %% use preexisting seed
  if (~exist(file512seed) & ~exist(fileNseed))
    disp(['You let oldseedflag true, but neither ' file512seed ' nor ' fileNseed ' exists.']);
    returnflag=true;
    return;
  end
  if exist(fileNseed)  %% use preexisting fileNseed
    if matlabflag
      load(fileNseed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', fileNseed, 'randamp', 'randphs');
    end
    clear fileNseed;
  elseif (Nmode_p < 512)  %% use preexisting file512seed for Nmode_p <= 512
    if matlabflag
      load(file512seed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', file512seed, 'randamp', 'randphs');
    end
    clear file512seed;  %% save memory 
    randamp = reshape(randamp, 512, 512, 257); %% into 3D
    randphs = reshape(randphs, 512, 512, 257); %% into 3D

    istart  = (512-Nmode_p)/2 + 1;  %% x & y starting index of 512*512*257 seed to use
    iend    = istart + Nmode_p -1;  %% x & y ending index of 512*512*257 seed to use
    %% Following two lines reduce the size of matrix.
    randamp = randamp(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
    randphs = randphs(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
    randamp = reshape(randamp, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
    randphs = reshape(randphs, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
  elseif (Nmode_p > 512)  %% use preexisting 512seed, and generate & attach higher-k seeds.
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
    randphs__ = unifrnd(0,2*pi, [Nmissing,1]);   %% new seed (phase)
    istart = (Nmode_p-512)/2 + 1;   %% x & y starting index to host 512*512*257 seed
    iend   = istart    + 512 - 1;   %% x & y ending index to host 512*512*257 seed
    izstart = Nc_p-257+1;           %% z starting index to host 512*512*257 seed

    %% host subgaussseed512.matbin seed
    randamp_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randamp, 512, 512, 257);
    randphs_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randphs, 512, 512, 257);
    clear randamp randphs;  %% save memory
    
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

    randamp = reshape(randamp_, Nmode_p*Nmode_p*Nc_p, 1);
    randphs = reshape(randphs_, Nmode_p*Nmode_p*Nc_p, 1);
    clear randamp_ randphs_;
  end
else  %% generate completely new seed
  if record512companionflag
    if (Nmode_p < 512)
      randamp = raylrnd(1,      [512*512*257,1]);   %% new seed (amplitude)
      randphs = unifrnd(0,2*pi, [512*512*257,1]);   %% new seed (phase)
      if matlabflag
        save(file512seed, 'randamp', 'randphs', '-v6');
      else
        save('-mat-binary', file512seed, 'randamp', 'randphs');
      end
      randamp = reshape(randamp, 512, 512, 257); %% into 3D
      randphs = reshape(randphs, 512, 512, 257); %% into 3D

      istart  = (512-Nmode_p)/2 + 1;  %% x & y starting index of 512*512*257 seed to use
      iend    = istart + Nmode_p -1;  %% x & y ending index of 512*512*257 seed to use
      %% Following two lines reduce the size of matrix.
      randamp = randamp(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
      randphs = randphs(istart:iend, istart:iend, 257-Nc_p+1:257);  %% final seed in 3D
      randamp = reshape(randamp, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
      randphs = reshape(randphs, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
    elseif (Nmode_p > 512)
      randamp_ = raylrnd(1,      Nmode_p*Nmode_p*Nc_p, 1);   %% new seed (amplitude)
      randphs_ = unifrnd(0,2*pi, Nmode_p*Nmode_p*Nc_p, 1);   %% new seed (phase)
      randamp_ = reshape(randamp_, Nmode_p, Nmode_p, Nc_p);  %% into 3D
      randphs_ = reshape(randphs_, Nmode_p, Nmode_p, Nc_p);  %% into 3D

      istart   = (Nmode_p-512)/2 + 1;  %% x & y starting index of randamp_ & randphs_
      iend     = istart + 512 -1;      %% x & y ending index of randamp_ & randphs_
      randamp  = randamp_(istart:iend, istart:iend, Nc_p-257+1:Nc_p);  %% extract 512seed
      randphs  = randphs_(istart:iend, istart:iend, Nc_p-257+1:Nc_p);  %% extract 512seed
      randamp  = reshape(randamp, 512*512*257, 1);  %% into 1D array
      randphs  = reshape(randphs, 512*512*257, 1);  %% into 1D array
      if matlabflag
        save(file512seed, 'randamp', 'randphs', '-v6');
      else
        save('-mat-binary', file512seed, 'randamp', 'randphs');
      end
      clear randamp randphs; %% save memory
      randamp = reshape(randamp_, Nmode_p*Nmode_p*Nc_p, 1);
      randphs = reshape(randphs_, Nmode_p*Nmode_p*Nc_p, 1);
      clear randamp_ randphs_; %% save memory
    end
  else
    randamp = raylrnd(1,      [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (amplitude)
    randphs = unifrnd(0,2*pi, [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (phase)
  end
end

if recordseedflag
  if matlabflag
    save(fileNseed, 'randamp', 'randphs', '-v6');
  else
    save('-mat-binary', fileNseed, 'randamp', 'randphs');
  end
end
