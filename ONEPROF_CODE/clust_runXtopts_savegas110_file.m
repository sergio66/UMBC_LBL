%% this simply does all wavenumbers for g1

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 1

%% file will contain AB CDEFG HI  which are gasID, wavenumber, temp offset   
%%                   12 34567 89
%% where gasID = 01 .. 99,   HI = 1 .. 11 (for Toff = -5 : +5) and wavenumber = 00050:99999

gasIDlist = load('g1_ir_list.txt');
XJOB = num2str(gasIDlist(JOB));
if length(XJOB) == 8
  XJOB = ['0' num2str(gasIDlist(JOB))];
end
  
Sgid     = str2num(XJOB(1:2));
Schunk   = str2num(XJOB(3:7));  
Stoffset = str2num(XJOB(8:9)); Spp = Stoffset - 6;

fprintf(1,'JOB String = %s    parsed to gid = %2i chunk = %5i Stoffset = %2i \n',XJOB,Sgid,Schunk,Stoffset);

nbox = 5;
pointsPerChunk = 10000;
gases = [1];

poffset = [0.1, 1.0, 3.3, 6.7, 10.0];

%% in /home/sergio/HITRAN2UMBCLBL      refproTRUE.mat -> refprof_usstd16Aug2010_lbl.mat
%% load /home/sergio/abscmp/refproTRUE.mat
load /home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat

addpath /home/sergio/SPECTRA
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

freq_boundaries_g110

cd /home/sergio/SPECTRA

%% qtips04.m : for water the isotopes are 161  181  171  162  182  172
%% so HDO = HOD = 162

%% don't need concept of JOB for G1 (so few chunks) but let us prototype anyway
fmin0 = fmin;

if Schunk >= fmin0
  fmin = Schunk;
else
  Schunk
  disp('the start wavnumber is SMALLER than fmain = 1105 cm-1 so quit')
  return
end

%% qtips04.m : for water the isotopes are 161  181  171  162  182  172
%% so HDO = HOD = 162
while fmin <= wn2
  fmax = fmin + dv;
  for pp = Spp
    gg = 1;
    gg = Sgid;
    fprintf(1,'gas freq = %3i %6i \n',gg,fmin);
    gasid = gg;  

    tprof = refpro.mtemp + pp*10;
    for mm = 1 : 5

      iYes = findlines_plot(fmin-25,fmax+25,1); 

      fout = [dirout '/stdH2OALL' num2str(fmin)];
      fout = [fout '_' num2str(gg) '_' num2str(pp+6) '_' num2str(mm) '.mat'];
      if exist(fout,'file') == 0 & iYes > 0
        toucher = ['!touch ' fout]; %% do this so other runs go to diff chunk 
        eval(toucher);
        profile = [(1:100)' refpro.mpres refpro.gpart(:,gg)*poffset(mm)  tprof refpro.gamnt(:,gg)]';
        fip = ['IPFILES/std_gx' num2str(gg) 'x_' num2str(pp+6) '_' num2str(mm)];
        fid = fopen(fip,'w');
        fprintf(fid,'%3i %10.8f %10.8f %7.3f %10.8e \n',profile);
        fclose(fid);

        topts.CKD = -1;
        topts.which_isotope = [1 2 3 4 5 6 7];; %% use isotopes all
        topts.which_isotope = [1 2 3   5 6 7];; %% use all EXCEPT isotope 4 = HDO
        topts.which_isotope = [-1 4];;          %% use all isotopes except 4 = HDO
        topts.which_isotope = [0   ];;          %% use all isotopes
        [w,d] = run8water(gasid,fmin,fmax,fip,topts);  

        %fout = [dirout '/stdH2OALL' num2str(fmin)];
        %fout = [fout '_' num2str(gg) '_' num2str(pp+6) '_' num2str(mm) '.mat'];
        saver = ['save ' fout ' w d '];
        eval(saver);
      elseif exist(fout,'file') > 0 & iYes > 0
        fprintf(1,'file %s already exists \n',fout);
      elseif exist(fout,'file') == 0 & iYes < 0
        fprintf(1,'no water lines for chunk starting %8.6f \n',fmin);
      end
    end             %% loop over partial pressure
  end               %% loop over temperature (1..11)
  fmin = fmin + dv;
  %% one chunk is enough
  return
end                 %% loop over freq

