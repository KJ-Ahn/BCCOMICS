%% Create/Use/Shuffle/Add random numbers for initial condition
%% 'subgaussseed512.matbin' is used as the seed for any Ncell_p<=512 cases,
%% when oldseedflag is true. For Ncell_p>512 cases, if oldseedflag is still true,
%% 'subgaussseed512.matbin' is used as part of the seed.

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
      load(file512seed, '-mat', 'randamp', 'randphs')
    else
      load('-mat-binary', file512seed, 'randamp', 'randphs')
    end
    istart = (512-Nmode_p)/2 + 1;
    iend   = istart + Nmode_p -1;
    randamp_ = randamp(istart:iend, istart:iend, 257-Nc_p+1:257);
    randphs_ = randphs(istart:iend, istart:iend, 257-Nc_p+1:257);
  elseif exist(fileNseed)  %% use preexisting seed, and also generate & attach higher-k seeds.
    if matlabflag
      load(fileNseed, '-mat', 'randamp', 'randphs')
    else
      load('-mat-binary', fileNseed, 'randamp', 'randphs')
    end
  elseif ~exist(fileNseed)  %% use preexisting seed, and also generate & attach higher-k seeds.
    if matlabflag
      load(file512seed, '-mat', 'randamp', 'randphs')
    else
      load('-mat-binary', file512seed, 'randamp', 'randphs')
    end
    Nmissing = Nmode_p*Nmode_p*Nc_p - 512*512*257;
    randamp_  = zeros(Nmode_p, Nmode_p, Nc_p);
    randphs_  = zeros(Nmode_p, Nmode_p, Nc_p);
    randamp__ = raylrnd(1,      [Nmissing,1]);
    randphs__ = unifrnd(0.2*pi, [Nmissing,1]);
    istart = (Nmode_p-512)/2 + 1;
    iend   = istart    + 512 - 1;
    izstart = Nc_p-257+1;

    Ntemp     = Nmode_p + Nmode_p + Nc_p-257;
    randamp_(1:Nmode_p, 1:Nmode_p, 1:Nc_p-257) = reshape(randamp__(1:Ntemp), Nmode_p, Nmode_p, Nc_p-257);
    Ntemp_old = Ntemp;
    Ntemp     = istart-1 + Nmode_p + 257;
    randamp_(1:istart-1, 1:Nmode_p, izstart:Nc_p) = reshape(randamp__(Ntemp_old+1:Ntemp_old+Ntemp), istart-1, Nmode_p, 257);
    randamp_(iend+1:Nmode_p, 1:Nmode_p, izstart:Nc_p) = reshape(randamp__(1:Ntemp), istart-1, Nmode_p, 257);
    iend+1:Nmode_p
    1:257-Nc_p
  end
else
    
    
    
end
