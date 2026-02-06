function [line,hitran_version,hlist_qtips] = ...
            hitread(start,stop,strengthM,gasID,filename,iCO2in);

%%    hlist_qtips tells you which HXX versions use the old qtips, and which
%%    use the new Lagrange stuff from the H04 database
% function line=hitreadNEW(start,stop,strengthM,gasID,filename,iCO2);
% this calls read_hitran, and adds on an extra field : LINCT which is
% the number of lines read in

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% warning .. need to agree with  %%%%
%%%%   qoldVSqnew_fast
%%%%   initializeABCDG_qtips
%%%% so do 
%%%%   grep -in 'length(intersect(hitran_version,hlist_qtips))' *.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simple version, w/o symbolic link checking/making is hitread_no_symbolic.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it differentiates and calls the right MEXED HITRAN readers
% >>>> when updating HITRAN versions, do a search for "UPDATE"

% addpath /asl/matlib/aslutil
% addpath /asl/matlib/read_hitr06/
% addpath /home/sergio/git/SPECTRA/read_hitr06/

addpath /home/sergio/git/sergio_matlib/matlib/aslutil
addpath /home/sergio/git/sergio_matlib/matlib/read_hitr06/
addpath /home/sergio/git/UMBC_LBL/read_hitr06/
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/read_hitr06/

current_dir = pwd;

blah = findstr('/',filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename
% HITRAN0 = '/asl/data/hitran/h16.by.gas/g2.dat';
% HITRAN1 = '/asl/rta/hitran/h16.by.gas//g2.dat';  %% this would totally fool it (ie the //) so WATCH OUT
% HITRAN1 = /umbc/xfs3/strow/asl/rta/hitran/h24.by.gas/g2.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% look at the fourth occurence or fifth occurence
if length(blah) == 5
  %disp('usual data version eg /asl/rta/hitran/h16.by.gas/g2.dat')
  %disp('usual data version eg /asl/data/hitran/h16.by.gas/g2.dat')
  disp('usual data version eg /asl/data/hitran/h20.by.gas/g2.dat')
  hitran_version = filename(blah(4)+1:blah(4)+3);  %% eg /asl/data/hitran/h16.by.gas/g2.dat
  iHiTemp = -1;
elseif length(blah) == 6
  %disp('HITEMP version eg /asl/data/hitran/HITEMP/h16.by.gas/g2.dat')
  disp('HITEMP version eg /asl/data/hitran/HITEMP/h20.by.gas/g2.dat')
  hitran_version = filename(blah(5)+1:blah(5)+3);  %% eg /asl/data/hitran/HITEMP/h16.by.gas/g2.dat
  iHiTemp = +1;
end

%% look at the fourth occurence or fifth occurence
if length(blah) == 8
  %disp('usual data version eg /asl/rta/hitran/h16.by.gas/g2.dat')
  %disp('usual data version eg /asl/data/hitran/h16.by.gas/g2.dat')
  disp('usual data version eg /umbc/xfs3/strow/asl/rta/hitran/h16.by.gas/g2.dat')
  hitran_version = filename(blah(7)+1:blah(7)+3);  %% eg /asl/data/hitran/h16.by.gas/g2.dat
  iHiTemp = -1;
elseif length(blah) == 9
  %disp('HITEMP version eg /asl/data/hitran/HITEMP/h16.by.gas/g2.dat')
  %disp('HITEMP version eg /asl/data/hitran/HITEMP/h20.by.gas/g2.dat')
  disp('HITEMP version eg /umbc/xfs3/strow/asl/rta/hitran/HITEMP/h20.by.gas/g2.dat')
  hitran_version = filename(blah(8)+1:blah(8)+3);  %% eg /asl/data/hitran/HITEMP/h16.by.gas/g2.dat
  iHiTemp = +1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UPDATE
old_reader = {'h92','h96','h98','h2k'};
new_reader = {'h04','h08','h12','h16','g15'};
new_reader = {'h04','h08','h12','h16','g15','h17','h18'};              %% h17 and h18 are junk so we can read in
                                                                       %% old/new HITRAN LM databases for CO2
new_reader = {'h04','h08','h12','h16','g15','h17','h18','h20','h24'};  %% h17 and h18 are junk so we can read in
                                                                       %% old/new HITRAN LM databases for CO2

iOld = -1;
iNew = -1;
%% first try REGULAR
for ii = 1 : length(old_reader)
  if strcmp(old_reader(ii),hitran_version)
    iOld = +1;
  end
end
for ii = 1 : length(new_reader)
  if strcmp(new_reader(ii),hitran_version)
    iNew = +1;
  end
end
if iOld == -1 & iNew == -1
  filename
  error('whoops : have NOT been able to figure what HITRAN version you want!')
end

if iOld > 0
  %fprintf(1,'--> using old MEX reader for HITRAN version %s : %s\n',hitran_version,filename)
  %line = read_hitran(start,stop,strengthM,gasID,filename);
  fprintf(1,'--> using old MATLAB reader for HITRAN version %s : %s \n',hitran_version,filename)
  line = read_hitranOLD_H92_H2k(start,stop,strengthM,gasID,filename);
elseif iOld < 0
  fprintf(1,'--> using new MEX reader for HITRAN version %s : %s\n',hitran_version,filename)
  % cd /asl/matlab2012/read_hitr06
  % cd /asl/matlab/read_hitr06
  %%%% addpath /asl/matlab2012/read_hitr06
  % cd /home/sergio/SPECTRA/read_hitr06
  %%%% addpath /home/sergio/SPECTRA/read_hitr06
  iFoastOrSlow = +1;  %% mex basd reader
  iFoastOrSlow = -1;  %% matlab reader
  if iFoastOrSlow > 0
    line = read_hitran(start,stop,strengthM,gasID,filename);
  else
    disp('WARNING : UMBC_LBL/hitread.m : using slow matlab reader for getting hitran data : read_hitran2')
    line = read_hitran2(start,stop,strengthM,gasID,filename);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line.linct=length(line.wnum);

line.iso=line.iso';
line.wnum=line.wnum';
line.stren=line.stren';
line.tprob=line.tprob';

line.abroad=line.abroad';
line.sbroad=line.sbroad';

line.els=line.els';
line.abcoef=line.abcoef';
line.tsp=line.tsp';

if iOld == -1 & iNew == +1
  line.iusgq=line.iusgq';
  line.ilsgq=line.ilsgq';
end

line.gasid=line.igas';

if line.linct > 0
  line.igas=line.gasid(1);
end

cder = ['cd ' current_dir];
eval(cder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the "mass.dat" symbolic link
%% make sure you point mass.dat --> massXY.dat if XY = 92,96,98,2k,04

iLinkMASS = -1;

% UPDATE
if hitran_version == 'h92'
  disp('oops : for mass.dat use H96 masses for H92');
  linker = '96';
elseif hitran_version == 'h96'
  linker = '96';
elseif hitran_version == 'h98'
  linker = '98';
elseif hitran_version == 'h2k'
  linker = '00';
  linker = '98';   %% yuk  
elseif hitran_version == 'h04'
  linker = '04';
elseif hitran_version == 'h08'
  linker = '08';
elseif hitran_version == 'h12'
  linker = '12';
elseif hitran_version == 'g15'
  disp('linking GEISA 2015 to H2016 mass,qtips etc')
  linker = '16';
elseif hitran_version == 'h16'
  linker = '16';
elseif hitran_version == 'h17'
  linker = '16';
elseif hitran_version == 'h18'
  linker = '16';
elseif hitran_version == 'h20'
  linker = '20';
elseif hitran_version == 'h24'
  linker = '24';
else
  error('hitread.m unknown hitran version for mass')
end

symlinker = ['!/bin/ln -s MASS_ISOTOPES/mass' linker '.dat mass.dat']; 

ee = exist('mass.dat','file');
if ee == 0
  fprintf(1,'mass.dat DNE .. making VERS %s symbolic link! \n',linker);
  %fprintf(1,'symbolic link command : %s \n',symlinker);
  iLinkMASS = +1;  
  eval(symlinker);
else
  cd /home/sergio/SPECTRA
  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  frand   = mktemp('ugh');
  lser = ['!ls -lt mass.dat >& ' frand];
  eval(lser)
  fid = fopen(frand,'r');
  ugher = fscanf(fid,'%s');
%  fprintf(1,'debugging hitread.m .... mass linker %s \n',ugher)
  fclose(fid);
  rmer = ['!/bin/rm ' frand];
  eval(rmer);
  eee = strcmp('lrwxrwxrwx',ugher(1:10));
  if eee == 0
    error('oh oh mass.dat is NOT a symbolic link; do not want to delete it!!!')
  else
    %%check to see if the symbolic link is for the right version
    tata = findstr('->mass',ugher);
    oldvers = ugher(tata+6:tata+7);
    fff = strcmp(oldvers,linker);
    %% if fff == 1, the old version == current wanted version ==> do nothing!
    if fff == 0
      boo = ['MASS_ISOTOPES/mass' linker '.dat'];
      eex = exist(boo,'file');
      if eex ~= 2
        fprintf(1,'trying to make symbolic link to %s but file DNE \n',boo)
        error('return in hitread.m');
      else
        fprintf(1,'----> updating VERS %s symbolic link for mass.dat! \n',linker)
        rmer = ['!/bin/rm mass.dat']; eval(rmer);
        eval(symlinker);
        iLinkMASS = +1;  
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the "qtips.m" symbolic link if HITRAN versions are H92,96,98
%% if HITRAN versions is H2k, symbolically link to H98
%% else use the new f77 qtips code for H04 onwards
%% make sure you point qtips.m --> qtipsXY.dat  if XY = 92,96
%%               point qtips.m --> qtips98.dat  if XY = 98,2k
%%               point qtips.m --> qtips04.dat  if XY = 04

iLinkQTIPS = -1;  

% UPDATE
if hitran_version == 'h92'
  linker = '92';
elseif hitran_version == 'h96'
  linker = '96';
elseif hitran_version == 'h98'
  linker = '98';
elseif hitran_version == 'h2k'
  linker = '04';
  linker = '00';
  linker = '98';
elseif hitran_version == 'h04'
  linker = '04';
  %linker = '98';   %% yuk  
elseif hitran_version == 'h08'
  linker = '08';
elseif hitran_version == 'h12'
  linker = '12';
elseif hitran_version == 'h16'
  linker = '16';
elseif hitran_version == 'h20'
  linker = '20';
elseif hitran_version == 'h24'
  linker = '24';
end

% UPDATE
all_hitran_versions = {'h92','h96','h98','h2k','h04','h08','h12','h16','h20','h24'};

symlinker = ['!/bin/ln -s qtips' linker '.m qtips.m']; 

ee = exist('qtips.m','file');  

if (ee == 0)   %% ee == 0 ==> does not exist
  if (length(intersect(hitran_version,all_hitran_versions)) == 1)
    symlinker = ['!/bin/ln -s qtips' linker '.m qtips.m']; 
    fprintf(1,'qtips.m  DNE .. making VERS %s symbolic link! \n',linker);
    boo = ['qtips' linker '.m'];
    eex = exist(boo,'file');
    if eex ~= 2
      fprintf(1,'trying to make symbolic link to %s but file DNE \n',boo)
      error('return in hitread.m');
    else
      eval(symlinker);
      iLinkQTIPS = +1;  
    end
  end
end

%%%%%%%%%%%%%%%%%% ------------->> << ----------------------
%% the Hitran versions that can use old polynom qtips is hardcoded here
%hlist_qtips = {'h92','h96','h98','h2k','h04'};    %% these can use old qtips
% UPDATE
hlist_qtips = {'h92','h96','h98','h2k'};                 %% these can use old qtips
hlist_qnew = {'h04','h08','h12','h16','h20','h24'};      %% these can use new lagrange

disp(' run8 will use old polynomial qtips for these HITRAN versions ');
for lll = 1 : length(hlist_qtips)
  fprintf(1,'  %s  ',hlist_qtips{lll})
end
fprintf(1,' \n')
disp(' run8 will use new lagrange qtips for these HITRAN versions ');
for lll = 1 : length(hlist_qnew)
  fprintf(1,'  %s  ',hlist_qnew{lll})
end
fprintf(1,' \n')
%%%%%%%%%%%%%%%%%% ------------->> << ----------------------

if (ee ~= 0) %% ee ~= 0 ==> exists; check to see if it is version we want
  if (length(intersect(hitran_version,all_hitran_versions)) == 1)
    cd /home/sergio/SPECTRA
    cd /home/sergio/SPECTRA    
    randstr= [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
    frand   = ['ugh' randstr];
    frand   = mktemp('ugh');
    lser = ['!ls -lt qtips.m >& ' frand];
    eval(lser)
    fid = fopen(frand,'r');
    ugher = fscanf(fid,'%s');
    fclose(fid);
    rmer = ['!/bin/rm ' frand];
    eval(rmer);
    eee = strcmp('lrwxrwxrwx',ugher(1:10));
%    fprintf(1,'debugging hitread.m .... qtips linker %s \n',ugher)    
    if eee == 0
      error('oops qtips.m is NOT a symbolic link; do not want to delete it!!!')
    else
      %%check to see if the symbolic link is for the right version
      tata = findstr('->qtips',ugher);
      oldvers = ugher(tata+7:tata+8);
      fff = strcmp(oldvers,linker);
      %% if fff == 1, the old version == current wanted version ==> do nothing!
      if fff == 0
        %% symbolic link is to a different version than current wanted version
        fprintf(1,'----> updating VERS %s symbolic link for qtips.m ! \n',linker)
        rmer = ['!/bin/rm qtips.m']; eval(rmer);
        %fprintf(1,'symbolic link command : %s \n',symlinker);
        eval(symlinker);
        iLinkQTIPS = +1;  
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hitran_version ~= 'h17' & hitra_version ~= 'h18'
  % UPDATE CO2_MATFILES
  update_co2_matfiles_hitread

  % UPDATE CO_MATFILES
  %update_co_matfiles_hitread  %% when I was trying CO line mixing
end
