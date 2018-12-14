%% Analyze & plot figures
%% 1. Does quick and easy analysis, and dumps itous output to 'IC_stats.dat'
%% 2. Plots some useful slice maps under directory /figs.
%%
%% Edit this to your taste for more/better analysis. Especially 'Ns' value
%% below is worth changing for visibility of quiver plots, and 'pkgdir'
%% need to be corrected.


pkgdir = '/home/kjahn/BCCOMICS'  %% Of course use your own package directory !!

%%%%% Option 1 -- begin
%%Consts_Conversions; %% need to copy here from src directory
%%params;
%%LCDM;  %% use Cosmology script of your own, or just define h (e.g. h=0.703).
%%%%% Option 1 -- end
%%%%% Option 2 -- begin
h=0.703;
%%%%% Option 2 -- end
params_patch;

if (exist('OCTAVE_VERSION','builtin'))
  matlabflag=false;
  disp('----------------run on OCTAVE----------------');
else
  matlabflag=true;
  disp('----------------run on MATLAB----------------');
end

if ~matlabflag
  setenv('GNUTERM','qt')  %% if the installed gnuplot is using QT instead of x11 as default. Comment this out and choose the other if error occurs when plotting.
%%  setenv('GNUTERM','x11')  %% if the installed gnuplot is using x11 instead of QT as default. Comment this out if error occurs when plotting.
end

%% If on terminal without x11 or qt capability, uncomment the following:
figure('visible','off');

%%%%%%%%%%%%%%%%%%%% simple analysis for sanity check %%%%%%%%%%%%%%%% begin
Npstr = num2str(Ncell_p);

fin = fopen('Units.txt','r');
fgets(fin);
fgets(fin);
fgets(fin);
Units = fscanf(fin, '%e %e %e');
fclose(fin);

zz = load('zz.dat');
%% redshift of initial condition
zred = zz(2);

fin = fopen('vc1','r');
vc1 = fread(fin, 'double');
vc1 = reshape(vc1, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvc1 = std(vc1(:))*Units(2);
avevc1 = sum(vc1(:))/Ncell_p^3*Units(2);
clear vc1;

fin = fopen('vc2','r');
vc2 = fread(fin, 'double');
vc2 = reshape(vc2, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvc2 = std(vc2(:))*Units(2);
avevc2 = sum(vc2(:))/Ncell_p^3*Units(2);
clear vc2;

fin = fopen('vc3','r');
vc3 = fread(fin, 'double');
vc3 = reshape(vc3, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvc3 = std(vc3(:))*Units(2);
avevc3 = sum(vc3(:))/Ncell_p^3*Units(2);
clear vc3;

fin = fopen('vb1','r');
vb1 = fread(fin, 'double');
vb1 = reshape(vb1, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvb1 = std(vb1(:))*Units(2);
avevb1 = sum(vb1(:))/Ncell_p^3*Units(2);
clear vb1;

fin = fopen('vb2','r');
vb2 = fread(fin, 'double');
vb2 = reshape(vb2, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvb2 = std(vb2(:))*Units(2);
avevb2 = sum(vb2(:))/Ncell_p^3*Units(2);
clear vb2;

fin = fopen('vb3','r');
vb3 = fread(fin, 'double');
vb3 = reshape(vb3, Ncell_p, Ncell_p, Ncell_p);
fclose(fin);

stdvb3 = std(vb3(:))*Units(2);
avevb3 = sum(vb3(:))/Ncell_p^3*Units(2);
clear vb3;

diary 'IC_stats.dat';
disp('-------------------------------');
disp(['IC at redshift ' num2str(zred)]);
disp('');
disp(['# of particles : ' Npstr '*' Npstr '*' Npstr]);
disp(['# of grid cells: ' Npstr '*' Npstr '*' Npstr]);
disp('');
disp(['unit of density          : ' num2str(Units(1)) ' g/cm^3']);
disp(['unit of peculiar velocity: ' num2str(Units(2)) ' cm/s']);
disp(['unit of specific energy  : ' num2str(Units(3)) ' cm^2/s^2']);
disp('');
disp('Average peculiar velocity of CDM (x,y,z) in km/s:');
disp([num2str(avevc1/1e5) ' ' num2str(avevc2/1e5) ' ' num2str(avevc3/1e5)]);
disp('Average peculiar velocity of baryon (x,y,z) in km/s:');
disp([num2str(avevb1/1e5) ' ' num2str(avevb2/1e5) ' ' num2str(avevb3/1e5)]);
disp('');
disp('standard deviation of peculiar velocity of CDM (x,y,z) in km/s:');
disp([num2str(stdvc1/1e5) ' ' num2str(stdvc2/1e5) ' ' num2str(stdvc3/1e5)]);
disp('standard deviation of peculiar velocity of baryon (x,y,z) in km/s:');
disp([num2str(stdvb1/1e5) ' ' num2str(stdvb2/1e5) ' ' num2str(stdvb3/1e5)]);
diary off;

%%%%%%%%%%%%%%%%%%%% simple analysis for sanity check %%%%%%%%%%%%%%%% end

%%%%%%%%%%%%%%%%%%%% useful plots; edit to your taste %%%%%%%%%%%%%%%% begin

if matlabflag
  load('4fig.matbin', '-mat', 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot') %%matlab
else
  load('-mat-binary', '4fig.matbin', 'xCDM_plane', 'xCDM_ex_plane', 'yCDM_plane', 'yCDM_ex_plane', 'Zc', 'Zb', 'Vc1', 'Vc2', 'Vc3', 'ZThc', 'Vb1', 'Vb2', 'Vb3', 'ZThb', 'Zeth', 'Ztemp', 'Zetot') %%octave
end

if matlabflag
  set(groot,'DefaultFigureColormap',jet);
end

%% extract cell indices, etc.
dirname = pwd;
Ndirname = length(dirname)
idx_L   = strfind(dirname,"Mpch_");
idx_i   = strfind(dirname,"_ic");
idx_j   = strfind(dirname,"_jc");
idx_k   = strfind(dirname,"_kc");
Nstr    = dirname(idx_L+4:idx_i-1);
icstr   = dirname(idx_i+3:idx_j-1 );
jcstr   = dirname(idx_j+3:idx_k-1 );
kcstr   = dirname(idx_k+3:Ndirname);

%% figure directory
figdir = 'figs'
mkdir(figdir);

%%%%%%%%%% Useful figures --------------------------- begin
p1    = imagesc(flipud(Zc)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'dc_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'dc_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

colormap(jet);
p1    = imagesc(flipud(Zb)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'db_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'db_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

p1    = imagesc(flipud(ZThc)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'dThc_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'dThc_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

p1    = imagesc(flipud(ZThb)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'dThb_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'dThb_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

p1   = imagesc(flipud(Zeth)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'Ethermbaryon_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'Ethermbaryon_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

p1    = imagesc(flipud(Ztemp)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'Temperature_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'Temperature_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

p1    = imagesc(flipud(Zetot)), colorbar;
set(gca, 'YTick', [], 'XTick', []);
colormap(jet);
pbaspect([1 1 1])
strout = [figdir '/' 'Etotbaryon_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'Etotbaryon_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

%% for padarray below.
if ~exist('padarray')
  addpath([pkgdir '/mfiles_for_matlab/ImagProc_Tool']);
end

%% colormap + quiver plot.
%% (would not work well on old versions(3.*.*) of Octave)
[X,Y]= meshgrid(1:Ncell_p+1);
p1=pcolor(X,Y,padarray(Zc,[1 1],'circular','post'));
shading flat;
axis nolabel;
colormap(jet);
pbaspect([1 1 1])
colorbar;
hold on;
%%%% if arrows become too densely packed, make it sparser by increasing Ns
Ns = 2; %% make quiver plot once every Nsparse elements
%%%% Vc1 is along the row, Vc2 is along the column, while
%%%% meshgrid defines X along the column and Y along the row.
%%%% So quiver(X,Y,Vc2,Vc1) is done..
p2 = quiver(X(1:Ns:Ncell_p,1:Ns:Ncell_p),Y(1:Ns:Ncell_p,1:Ns:Ncell_p),Vc2(1:Ns:Ncell_p,1:Ns:Ncell_p),Vc1(1:Ns:Ncell_p,1:Ns:Ncell_p),'filled');
set(p2,'color','black');
strout = [figdir '/' 'dc_vc_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'dc_vc_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

%% colormap + quiver plot.
%% (would not work well on old versions(3.*.*) of Octave)
p1=pcolor(X,Y,padarray(Zb,[1 1],'circular','post'));
shading flat;
axis nolabel;
colormap(jet);
pbaspect([1 1 1])
colorbar;
hold on;
p2 = quiver(X(1:Ns:Ncell_p,1:Ns:Ncell_p),Y(1:Ns:Ncell_p,1:Ns:Ncell_p),Vb2(1:Ns:Ncell_p,1:Ns:Ncell_p),Vb1(1:Ns:Ncell_p,1:Ns:Ncell_p),'filled');
set(p2,'color','black');
strout = [figdir '/' 'db_vb_i' icstr '_j' jcstr '_k' kcstr '.eps'];
print('-depsc', strout)
strout = [figdir '/' 'db_vb_i' icstr '_j' jcstr '_k' kcstr '.png'];
print('-dpng', strout)
close

