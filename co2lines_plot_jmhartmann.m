%%% this loads in the data from JM Hartmann IR 
%%  ~/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Data/BandInfo.dat
%% and loops thru to identify things by reverse-engineering co2lines_plot

hartmann = load('~/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Data/BandInfo.dat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% see HITRAN bd_vibs.for     for the spectroscopic notation  %%%%%%
%%%%%%            /asl/data/hitran/HITRAN2k/BD_Vibs.for           %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% this file finds the lines you want for the band in question
%%% copied from co2lines.m

gasID = 2;
%fnamePRE ='/asl/data/hitran/h92.by.gas/g';
fnamePRE ='/salsify/scratch4/h96.by.gas/g';
fnamePRE ='/asl/data/hitran/h98.by.gas/g';
fnamePRE = [hitranpath  hitran_version '.by.gas/g'];

fnamePOST = '.dat';
fnameIN = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];

start = 400;
stop = 2800;

%%line = hitread(start,stop,1.0e-28,gasID,hitlin_fname);
line = hitread(start,stop,0,gasID,hitlin_fname);

ggall   = [];
bandall = [];
[m,n] = size(hartmann);
str = ...
  '#  Iso    US     LS    Band    Numlines(Band)   NumLines(SoFar)   AllLines';
disp(str);
disp('-----------------------------------------------------------------------')
for ii = 1 : m
  xx = hartmann(ii,:);
  band = reverse_engineer_co2band(xx);
  bandall = [bandall band];
  %% now remove these lines
  gg = find(line.iso == xx(1) & line.iusgq == xx(2) &  line.ilsgq == xx(3));
  ggall = [ggall gg];
  semilogy(line.wnum,line.stren,'b.',line.wnum(ggall),line.stren(ggall),'r.')
  data = ...
    [ii xx(1) xx(2) xx(3) band length(gg) length(ggall) length(line.wnum)];
  fprintf(1,'%3i  %3i    %3i    %3i     %3i    %4i    %7i    %7i\n',data')
  pause;
  end
