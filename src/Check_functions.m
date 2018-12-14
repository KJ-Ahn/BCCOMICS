%% Some old versions of gnu octave has buggy ifftshift routine, so for
%% octave version older than 4.0.1, just use working one under the provided
%% directory. In case statistics package (for raylrnd) is not installed,
%% use provided statistics package.
%%
%% Also, interp2 is used for interpolating the transfer function, and 
%% some old version of octave interpn either requires meshgrid data for X & Y 
%% in interp2(X,Y,Z,x,y,interp2opt), or 'pchip' is not allowed as interp2opt, 
%% or even some old octave does not have interp2 function. 
%% Test this and if this is the case use the working interp2.m provides
%% in mfiles_for_octave.
%%
%% interpn is used for 3D interpolation when particlevelocity_accuracyflag
%% is true. Again, test this ('pchip' is not implemented yet, so use 'linear'
%% or 'spline' for interpolation method) and if error occurs let bccomics
%% use working interpn.m. For safety 'linear' is preferred: not smooth but
%% does not have some unlucky weird behavior of 'spline' method, and still
%% this is better than 1st order LPT that uses Eulerian time derivative for
%% velocity.
%%
%% All this can be avoided by upgrading to most recent
%% octave version and installing octave-statistics package.
if ~matlabflag
  if compare_versions(OCTAVE_VERSION,'4.0.1','<')
    %% Messages "warning: function * shadows ..." should be welcomed.
    addpath([pkgdir '/mfiles_for_octave']); 
  end
  try  %% test system interp2
    a=[1:3]';
    b=1:2;
    c=a*b;
    d=interp2(b,a,c,1.5,1.5,interp2opt);
  catch  %% when error occurs in interp2 use working interp2 instead 
    addpath([pkgdir '/mfiles_for_octave']); 
  end
  try  %% test system interpn
    x = -1:1;
    y = x;
    z = x;
    f = @(x,y,z) x.^2 - y - z.^2;
    [xx, yy, zz] = meshgrid (x, y, z);
    v = f (xx,yy,zz);
    xi = -1:0.1:1;
    yi = xi;
    zi = xi;
    [xxi, yyi, zzi] = ndgrid (xi, yi, zi);
    vi = interpn (x, y, z, v, xxi, yyi, zzi, interpnopt);
  catch  %% when error occurs in interpn use working interpn instead 
    addpath([pkgdir '/mfiles_for_octave']); 
  end
  try  %% test if padarray exists
    a=[1 2;3 4];
    aa=padarray(a,[1 1],'circular','post');
  catch  %% when error occurs guide for installation
    addpath([pkgdir '/mfiles_for_octave']); 
  end
  if ~exist('raylrnd')
    addpath([pkgdir '/statistics-1.3.0/inst']);
  end
else
  %% mfiles_for_matlab provides GPL licensed scripts that do what
  %% MATLAB Statistics Toolbox and Image Processing Toolbox would do
  %% for bccomics_setup and bccomics. Still, you may like to install
  %% these toolboxes for efficiency.
  if (~exist('raylrnd') || ~exist('unifrnd'))
    %%disp('Statistics Toolbox need be installed. Rerun after installation.');
    %%returnflag = true;
    %%return;
    addpath([pkgdir '/mfiles_for_matlab/Stat_Tool']);
  end
  if ~exist('padarray')
    %%disp('Image Processing Toolbox need be installed. Rerun after installation.');
    %%returnflag = true;
    %%return;
    addpath([pkgdir '/mfiles_for_matlab/ImagProc_Tool']);
  end
end
