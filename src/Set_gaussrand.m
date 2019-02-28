%% Create/Use/Shuffle/Add random numbers for initial condition.
%% First search for fileNseed and use it.
%% In case fileNseed doesnt exist but file512seed does,
%% 'subgaussseed512.matbin' is used as the seed for any Ncell_p<=512 cases,
%% when oldseedflag is true. For Ncell_p>512 cases, if oldseedflag is still true,
%% 'subgaussseed512.matbin' is used as part of the seed and missing seeds are
%% generated and attached.
%% Seeds are arranged such that resolution studies can be done.
%%
%% The end result is randamp and randphs, both as 1D arrays of Nmode_p*Nmode_p*Nc_p.

Nold  = Noldseed;
Ncold = Nold/2+1;

if (mod(Nold,4) ~= 0)
    disp('mod(Noldseed,4) should be zero. Check params_patch.m');
    return;
    returnflag = true;
end

if oldseedflag  %% use preexisting seed
  fileNoldseed = [diroldseed '/subgaussseed' num2str(Noldseed) '.matbin'];
  fileNseed    = [diroldseed '/subgaussseed' num2str(Nmode_p) '.matbin'];
  if (~exist(fileNoldseed) & ~exist(fileNseed))
    disp(['You let oldseedflag true, but neither ' fileNoldseed ' nor ' fileNseed ' exists.']);
    disp(['Set oldseedflag to false in params_patch.m and rerun bccomics (will generate new seed)']);
    returnflag=true;
    return;
  end
  if exist(fileNseed)  %% use preexisting fileNseed
    if matlabflag
      load(fileNseed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', fileNseed, 'randamp', 'randphs');
    end
    clear fileNseed;  %% save memory 
  elseif (Nmode_p < Nold)  %% use part of preexisting fileNoldseed for Nmode_p < Nold
    if matlabflag
      load(fileNoldseed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', fileNoldseed, 'randamp', 'randphs');
    end
    clear fileNoldseed;  %% save memory 
    randamp = reshape(randamp, Nold, Nold, Ncold); %% into 3D
    randphs = reshape(randphs, Nold, Nold, Ncold); %% into 3D

    istart  = (Nold-Nmode_p)/2 + 1;  %% x & y starting index of Nold*Nold*Ncold seed to use
    iend    = istart + Nmode_p -1;  %% x & y ending index of Nold*Nold*Ncold seed to use
    %% Following two lines reduce the size of matrix.
    randamp = randamp(istart:iend, istart:iend, Ncold-Nc_p+1:Ncold);  %% final seed in 3D
    randphs = randphs(istart:iend, istart:iend, Ncold-Nc_p+1:Ncold);  %% final seed in 3D
    randamp = reshape(randamp, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
    randphs = reshape(randphs, Nmode_p*Nmode_p*Nc_p, 1);  %% into 1D array
  elseif (Nmode_p > Nold)  %% use preexisting Noldseed, and generate & attach higher-k seeds.
    if matlabflag
      load(fileNoldseed, '-mat', 'randamp', 'randphs');
    else
      load('-mat-binary', fileNoldseed, 'randamp', 'randphs');
    end
    clear fileNoldseed;
    Nmissing = Nmode_p*Nmode_p*Nc_p - Nold*Nold*Ncold; %% # of missing modes
    randamp_  = zeros(Nmode_p, Nmode_p, Nc_p);   %% seed to fill in
    randphs_  = zeros(Nmode_p, Nmode_p, Nc_p);   %% seed to fill in
    randamp__ = raylrnd(1,      [Nmissing,1]);   %% new seed (amplitude)
    randphs__ = unifrnd(0,2*pi, [Nmissing,1]);   %% new seed (phase)
    istart = (Nmode_p-Nold)/2 + 1;   %% x & y starting index to host Nold*Nold*Ncold seed
    iend   = istart    + Nold - 1;   %% x & y ending index to host Nold*Nold*Ncold seed
    izstart = Nc_p-Ncold+1;           %% z starting index to host Nold*Nold*Ncold seed

    %% host subgaussseedNold.matbin seed
    randamp_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randamp, Nold, Nold, Ncold);
    randphs_(istart:iend, istart:iend, izstart:Nc_p)  = reshape(randphs, Nold, Nold, Ncold);
    clear randamp randphs;  %% save memory
    
    %% attach - 1 (bottom)
    is    = 1;                           %% starting index in 1D seed array
    Ntemp = Nmode_p*Nmode_p*(Nc_p-Ncold);  %% # of modes being filled in
    ie    = is + Ntemp-1;                %% ending index in 1D seed array  
    randamp_(1:Nmode_p, 1:Nmode_p, 1:Nc_p-Ncold)        = reshape(randamp__(is:ie), Nmode_p, Nmode_p, Nc_p-Ncold);  %% fill in part of randamp_
    randphs_(1:Nmode_p, 1:Nmode_p, 1:Nc_p-Ncold)        = reshape(randphs__(is:ie), Nmode_p, Nmode_p, Nc_p-Ncold);  %% fill in part of randphs_

    %% attach - 2 (top left)
    is    = ie + 1;                  %% starting index in 1D seed array
    Ntemp = (istart-1)*Nmode_p*Ncold;  %% # of modes being filled in     
    ie    = is + Ntemp-1;            %% ending index in 1D seed array  
    randamp_(1:istart-1, 1:Nmode_p, izstart:Nc_p)     = reshape(randamp__(is:ie), istart-1, Nmode_p, Ncold);  %% fill in part of randamp_
    randphs_(1:istart-1, 1:Nmode_p, izstart:Nc_p)     = reshape(randphs__(is:ie), istart-1, Nmode_p, Ncold);  %% fill in part of randphs_

    %% attach - 3 (top right)
    is    = ie + 1;                     %% starting index in 1D seed array
    Ntemp = (Nmode_p-iend)*Nmode_p*Ncold; %% # of modes being filled in     
    ie    = is + Ntemp-1;               %% ending index in 1D seed array  
    randamp_(iend+1:Nmode_p, 1:Nmode_p, izstart:Nc_p) = reshape(randamp__(is:ie), Nmode_p-iend, Nmode_p, Ncold);  %% fill in part of randamp_
    randphs_(iend+1:Nmode_p, 1:Nmode_p, izstart:Nc_p) = reshape(randphs__(is:ie), Nmode_p-iend, Nmode_p, Ncold);  %% fill in part of randphs_

    %% attach - 4 (top middle front)
    is    = ie + 1;                         %% starting index in 1D seed array
    Ntemp = (iend-istart+1)*(istart-1)*Ncold; %% # of modes being filled in     
    ie    = is + Ntemp-1;                   %% ending index in 1D seed array  
    randamp_(istart:iend, 1:istart-1, izstart:Nc_p)   = reshape(randamp__(is:ie), iend-istart+1, istart-1, Ncold);  %% fill in part of randamp_
    randphs_(istart:iend, 1:istart-1, izstart:Nc_p)   = reshape(randphs__(is:ie), iend-istart+1, istart-1, Ncold);  %% fill in part of randphs_

    %% attach - 5 (top middle back)
    is    = ie + 1;                             %% starting index in 1D seed array
    Ntemp = (iend-istart+1)*(Nmode_p-iend)*Ncold; %% # of modes being filled in     
    ie    = is + Ntemp-1;                       %% ending index in 1D seed array  
    randamp_(istart:iend, iend+1:Nmode_p, izstart:Nc_p)   = reshape(randamp__(is:ie), iend-istart+1, Nmode_p-iend, Ncold);  %% fill in part of randamp_
    randphs_(istart:iend, iend+1:Nmode_p, izstart:Nc_p)   = reshape(randphs__(is:ie), iend-istart+1, Nmode_p-iend, Ncold);  %% fill in part of randphs_

    randamp = reshape(randamp_, Nmode_p*Nmode_p*Nc_p, 1);
    randphs = reshape(randphs_, Nmode_p*Nmode_p*Nc_p, 1);
    clear randamp_ randphs_;
  end
else  %% generate completely new seed
  randamp = raylrnd(1,      [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (amplitude)
  randphs = unifrnd(0,2*pi, [Nmode_p*Nmode_p*Nc_p,1]);   %% new seed (phase)
end
