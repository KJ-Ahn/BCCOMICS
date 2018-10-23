## ratio of scale factors. Will be used for halo finding, etc.

tHi_a      = load('tHi_a.dat');
tHi_al     = load('tHi_alocal.dat');
tHi_aln    = load('tHi_alocal_n.dat');
tHi_aln123 = load('tHi_alocal_n123.dat');
tHi_aln4   = load('tHi_alocal_n4.dat');

NtHi = length(tHi_a(:,1));
tHi  = tHi_a(:,1);
a    = tHi_a(:,2);

al       = tHi_al(:,2:13);
aln5678  = tHi_aln(:,6:9);
aln123   = tHi_aln123(:,2:4);
aln4     = tHi_aln4(:,2);

Nn123 = length(aln123(:,1));
Nn4   = length(aln4(:,1));

ratio      = (a*ones(1,12))         ./al;
ration5678 = (a*ones(1,4))          ./aln5678;
ration123  = (a(1:Nn123)*ones(1,3)) ./aln123;
ration4    = a(1:Nn4)               ./aln4;

fout=fopen('a2alocal.dat','w');
dat = [tHi ratio];
fprintf(fout,'%e %e %e %e %e %e %e %e %e %e %e %e %e\n', dat');
fclose(fout);

fout=fopen('a2alocaln5678.dat','w');
dat_ = [tHi ration5678];
fprintf(fout,'%e %e %e %e %e\n', dat_');
fclose(fout);

fout=fopen('a2alocaln123.dat','w');
dat__ = [tHi(1:Nn123) ration123];
fprintf(fout,'%e %e %e %e\n', dat__');
fclose(fout);

fout=fopen('a2alocaln4.dat','w');
dat___ = [tHi(1:Nn4) ration4];
fprintf(fout,'%e %e\n', dat___');
fclose(fout);
