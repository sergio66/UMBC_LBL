%% this simply does all wavenumbers for g1

addpath /home/sergio/SPECTRA
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 1;
end

nbox = 5;
pointsPerChunk = 10000;
gases = [1];

poffset = [0.1, 1.0, 3.3, 6.7, 10.0];  %% we use mm === 2 below

gid = 1;

%% COMMENT
%% before we first called "rtp_prof_to_oneprof0.m" then it saved "refpro" into "oneprof.mat", which we loaded
%%  load /home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat         %% symbolic link
%%  load oneprof.mat
%% the main driver file was individual_clust_runXtopts_savegasN_file.m

%% now we call "rtp_prof_to_oneprof1.m which passes SHTUFF into "
%% [rtpProf] = rtp_prof_to_oneprof1(fname,5);
%%
%% rtpProf has ALREADY been made
load oneprof.mat
glist = rtpProf.glist;

freq_boundaries

cd /home/sergio/SPECTRA

%% qtips04.m : for water the isotopes are 161  181  171  162  182  172
%% so HDO = HOD = 162

%% don't need concept of JOB for G1 (so few chunks) but let us prototype anyway
fmin0 = fmin;

if JOB >= fmin0 & JOB <= fmax
  fmin = JOB;
else
  disp('huh')
  return
end

gasind = find(glist == gid);
if length(gasind) == 0
  glist
  gid
  error('cannot find gas in oneprof.mat')
end

%% qtips04.m : for water the isotopes are 161  181  171  162  182  172
while fmin <= wn2
  fmax = fmin + dv;
  for pp = 0
    gg = 1;
    fprintf(1,'gas freq = %3i %6i \n',gg,fmin);
    gasid = gg;  

    tprof = rtpProf.mtemp + pp*10;
    for mm = 2

      iYes = findlines_plot(fmin-25,fmax+25,1); 

      fout = [dirout '/profH2O' num2str(fmin)];
      fout = [fout '_' num2str(gg) '_' num2str(pp+6) '_' num2str(mm) '.mat'];
      if exist(fout,'file') == 0 & iYes > 0
        toucher = ['!touch ' fout]; %% do this so other runs go to diff chunk 
        eval(toucher);
        profile = [(1:100)' rtpProf.mpres rtpProf.gpart(:,gg)*poffset(mm)  tprof rtpProf.gamnt(:,gg)]';
        fip = ['IPFILES/std_gx' num2str(gg) 'x_' num2str(pp+6) '_' num2str(mm)];
        fid = fopen(fip,'w');
        fprintf(fid,'%3i %10.8f %10.8f %7.3f %10.8e \n',profile);
        fclose(fid);

        topts.CKD = -1;
%        topts.which_isotope = [1 2 3 4 5 6];; %% use isotopes all
%        topts.which_isotope = [1 2 3   5 6];; %% use all EXCEPT isotope 4 = HDO
%        topts.which_isotope = [-1 4];;        %% use all isotopes except 4 = HDO
        [w,d] = run8water(gasid,fmin,fmax,fip,topts);  

        saver = ['save ' fout ' w d profile '];
        eval(saver);
      elseif exist(fout,'file') > 0 & iYes > 0
        fprintf(1,'file %s already exists \n',fout);
      elseif exist(fout,'file') == 0 & iYes < 0
        fprintf(1,'no water lines for chunk starting %8.6f \n',fmin);
      end
    end             %% loop over partial pressure      mm
  end               %% loop over temperature (1..11)   pp
  fmin = fmin + dv;
  %% one chunk is enough
  return
end                 %% loop over freq

