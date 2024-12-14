%% this simply does all wavenumbers for xsec

% local running to test
% clustcmd -L clust_runXtopts_savexsecN_file.m '510173001'
%
% otherwise when happy
% clustcmd -q short -n 32 clust_runXtopts_savexsecN_file.m gN_ir_xseclist.txt
%
%% file will contain AB CDEFG HI  which are gasID, wavenumber, temp offset   
%%                   12 34567 89
%% where gasID = 01 .. 99,   HI = 1 .. 11 (for Toff = -5 : +5) and wavenumber = 00050:99999

Sgid     = str2num(JOB(1:2));
Schunk   = str2num(JOB(3:7));  
Stoffset = str2num(JOB(8:9)); Stt = Stoffset - 6;
  if Stt ~= 0
    error('Need Stt == 0')
  end

fprintf(1,'JOB String = %s    parsed to gid = %2i chunk = %5i Stoffset = %2i \n',JOB,Sgid,Schunk,Stoffset);

%% does the cross sections
%% see     head /asl/data/hitran/xsec98.ok/*.xsc | more

nbox = 5;
pointsPerChunk = 10000;

gases = 51:63; %%%%%%% <<<<<<<< need to change this as needed!
gases = Sgid;
gid   = Sgid;

%% load /home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat         %% symbolic link
load oneprof.mat

addpath /home/sergio/HITRAN2UMBCLBL/READ_XSEC/
addpath /home/sergio/SPECTRA
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

freq_boundaries

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

while fmin <= wn2
  fmax = fmin + dv;
  for tt = Stt
    wfreq0 = [];
    dkcomp0 = zeros(100,10000);
    for gg = 1 : length(gases)
      gasid = gases(gg);  

      gq = find(refpro.glist == gasid);
      tprof = refpro.mtemp + tt*10;
      profile = ...
        [(1:100)' refpro.mpres refpro.gpart(:,gq) tprof refpro.gamnt(:,gq)]';
      cd /home/sergio/SPECTRA

      [iYes,gf] = findxsec_plot(fmin,fmax,gasid); 
      fprintf(1,'gas freq = %3i %6i %3i %3i\n',gasid,fmin,iYes,tt);
      fout = [dirout '/std' num2str(fmin) '_' num2str(refpro.glist(gq)) '_' num2str(tt+6) '.mat'];

      ee = exist(fout,'file');
      if ee > 0
        fprintf(1,'%s already exists \n',fout);
      end

      %% see abscmp/xsectab25.m
      gamnt = profile(5,:);
      gamnt2d = ones(1e4,1) * gamnt;
      tp = profile(4,:);
      pL = profile(2,:);
      d = zeros(100,10000);
      if iYes > 0 & ee == 0
        toucher = ['!touch ' fout]; %% do this so other runs go to diff chunk  
        eval(toucher); 
        %[w,d] = calc_xsec(gf,fmin,fmax-dv,dv,tp,pL,1);  
        v1 = fmin; v2 = fmax-dv;
        dvv = dv/pointsPerChunk;
        nvpts = 1 + round((v2-v1)/dvv);
        %[gasid fmin fmax dv dvv nvpts]
        %[d,w] = calc_xsec(gasid,fmin,fmax-dvv,dvv,tp,pL,1);  
        [d,w] = calc_xsec(gasid,fmin,fmax-dvv,dvv,tp,pL);  
        d = d.*gamnt2d;
        d = d';    %%% need same dimensions as rest of gases!
      end

      if iYes > 0 & ee == 0
        saver = ['save ' fout ' w d '];
        eval(saver);
      end

    end               %% loop over gas
  end                 %% loop over temperature (1..11)
  fmin = fmin + dv;
  %% one chunk is enough
  return
end                   %% loop over freq

cd /home/sergio/HITRAN2UMBCLBL/MAKEIR_ALL_H12
