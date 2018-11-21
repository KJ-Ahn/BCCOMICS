%% Generates global redshifts. No special spacing rule needed, just sample
%% the way you like.
%% remove %% from below for your choice and run this script

%% piecewise linear
zglobal = [linspace(200,120,5)';  linspace(100,40,7)'; linspace(30,22,5)'; linspace(20,11,10)'; linspace(10,3,15)'];

%% fully linear sampling 
%% zglobal = linspace(200,5,40)';

%% fully log-linear sampling
%% zglobal = (10.^linspace(log10(200), log10(5), 50))';

fout = fopen('zglobal.dat','w');
fprintf(fout, '%.3f\n', zglobal');
fclose(fout);
