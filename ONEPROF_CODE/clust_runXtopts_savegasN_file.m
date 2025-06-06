%% this simply does all wavenumbers for gN 

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 1;
  JOB = 34;
end
%JOB = 1

%% file will contain AB CDEFG HI  which are gasID, wavenumber, temp offset   
%%                   12 34567 89
%% where gasID = 01 .. 99,   HI = 1 .. 11 (for Toff = -5 : +5) and wavenumber = 00050:99999

gasIDlist = load('gN_ir_list.txt');
XJOB = num2str(gasIDlist(JOB));
if length(XJOB) == 8
  XJOB = ['0' num2str(gasIDlist(JOB))];
end

Sgid     = str2num(XJOB(1:2));
Schunk   = str2num(XJOB(3:7));  
Stoffset = str2num(XJOB(8:9)); Stt = Stoffset - 6;   
  if Stt ~= 0
    error('Need Stt == 0')
  end

fprintf(1,'JOB String = %5i    parsed to gid = %2i chunk = %5i Stoffset = %2i \n',JOB,Sgid,Schunk,Stoffset);

nbox = 5;
pointsPerChunk = 10000;
gases = Sgid;

%% in /home/sergio/HITRAN2UMBCLBL      refproTRUE.mat -> refprof_usstd16Aug2010_lbl.mat
%% load /home/sergio/abscmp/refproTRUE.mat
%% load /home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat

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

addpath /home/sergio/SPECTRA
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

gg    = Sgid;
gasid = Sgid;  
gid   = Sgid;
freq_boundaries
wn1 = Schunk;
wn2 = Schunk + 25 - dv;

cd /home/sergio/SPECTRA

%% don't need concept of JOB for G1 (so few chunks) but let us prototype anyway
fmin0 = fmin;

if Schunk >= fmin0
  fmin = Schunk;
else
  Schunk
  disp('the start wavnumber is SMALLER than fmain = 1105 cm-1 so quit')
  return
end

gasind = find(glist == gg);
if length(gasind) == 0
  glist
  gg
  error('cannot find gas in oneprof.mat')
end

while fmin <= wn2
  fmax = fmin + dv;

  fprintf(1,'gas freq = %3i %6i \n',gg,fmin);

  for tt = Stt
    %tprof = refpro.mtemp + pp*10;
    tprof = rtpProf.mtemp + tt*10;

    iYes = findlines_plot(fmin-dv,fmax+dv,1); 

    fout = [dirout '/prof' num2str(fmin)];
    fout = [fout '_' num2str(gg) '_' num2str(tt+6) '.mat'];
    if exist(fout,'file') == 0 & iYes > 0
      toucher = ['!touch ' fout]; %% do this so other runs go to diff chunk 
      eval(toucher);
      %profile = [(1:100)' refpro.mpres refpro.gpart(:,gg)*poffset(mm)  tprof refpro.gamnt(:,gg)]';
      profile = [(1:length(rtpProf.mpres))' rtpProf.mpres rtpProf.gpart(:,gasind) tprof rtpProf.gamnt(:,gasind)]';
      fip = ['/home/sergio/SPECTRA/IPFILES/std_gx' num2str(gg) 'x_' num2str(tt+6) '.txt'];
      fid = fopen(fip,'w');
      fprintf(fid,'%3i %10.8f %10.8f %7.3f %10.8e \n',profile);
      fclose(fid);

      if gasid ~= 2
        [w,d] = run8(gasid,fmin,fmax,fip,topts);  
      else
        [w,d] = run8co2_linemixUMBC(gasid,fmin,fmax,fip,topts);  
      end

      saver = ['save ' fout ' w d profile '];
      eval(saver);
    elseif exist(fout,'file') > 0 & iYes > 0
      fprintf(1,'file %s already exists \n',fout);
    elseif exist(fout,'file') == 0 & iYes < 0
      fprintf(1,'no lines for chunk starting %8.6f \n',fmin);
    end
  end               %% loop over temperature (1..11)
  fmin = fmin + dv;
  %% one chunk is enough
  return
end                 %% loop over freq

