ipfile = 'IPFILES/n2one';

%[wx,dx] = run8(22,1305,2905,ipfile);      %% new lbl

f1 = 2300; f2 = 2400;
f1 = 2320; f2 = 2370;
f1 = 2130; f2 = 2530;

topts.O2O3N2continuumVers = 4;
[wx,dx4] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 4')

topts.O2O3N2continuumVers = 1;
[wx,dx1] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 1')

topts.O2O3N2continuumVers = 3;
[wx,dx3] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 3')

