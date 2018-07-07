function [line,hitran_version] = hitread(start,stop,strengthM,gasID,filename);

% function line=hitreadNEW(start,stop,strengthM,gasID);
% this calls read_hitran, and adds on an extra field : LINCT which is
% the number of lines read in

% it differentiates and calls the right MEXED HITRAN readers

current_dir = pwd;

blah = findstr('/',filename);
%%look at the fourth occurence
hitran_version = filename(blah(4)+1:blah(4)+3);

old_reader = {'h92','h96','h98','h2k'};
new_reader = {'h04','h08'};

iOld = -1;
iNew = -1;
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
  error('whoops : have NOT been able to figure what HITRAN version you want!')
  end

if iOld > 0
  fprintf(1,'using old MEX reader for HITRAN version %s \n',hitran_version)
  cd /asl/matlib/aslutil/
  line = read_hitran(start,stop,strengthM,gasID,filename);
else
  fprintf(1,'using new MEX reader for HITRAN version %s \n',hitran_version)
  cd /asl/matlib/read_hitr06
  line = read_hitran(start,stop,strengthM,gasID,filename);
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

line.iusgq=line.iusgq';
line.ilsgq=line.ilsgq';

line.gasid=line.igas';

if line.linct > 0
  line.igas=line.gasid(1);
  end

cder = ['cd ' current_dir];
eval(cder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the "mass.dat" symbolic link
%% make sure you point mass.dat --> massXY.dat if XY = 92,96,98,2k,04

if hitran_version == 'h92'
  disp('oops : for mass.dat use H96 masses for H92');
  linker = '96';
elseif hitran_version == 'h96'
  linker = '96';
elseif hitran_version == 'h98'
  linker = '98';
elseif hitran_version == 'h2k'
  linker = '00';
elseif hitran_version == 'h04'
  linker = '04';
elseif hitran_version == 'h08'
  linker = '08';
  end

symlinker = ['!/bin/ln -s mass' linker '.dat mass.dat']; 

ee = exist('mass.dat','file');
if ee == 0
  fprintf(1,'mass.dat not found .. making VERS %s symbolic link! \n',linker);
  %fprintf(1,'symbolic link command : %s \n',symlinker);
  eval(symlinker);
else
  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  lser = ['!ls -lt mass.dat >& ' frand];
  eval(lser)
  fid = fopen(frand,'r');
  ugher = fscanf(fid,'%s');
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
       %% symbolic link is to a different version than current wanted version
      fprintf(1,'updating VERS %s symbolic link for mass.dat! \n',linker)
      rmer = ['!/bin/rm mass.dat']; eval(rmer);
      eval(symlinker);
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the "qtips.m" symbolic link if HITRAN versions are H92,96,98
%% else we use the new f77 qtips code for H2k and H04
%% make sure you point qtips.m --> qtipsXY.dat  if XY = 92,96,98
%%               point qtips.m --> qtips04.dat  if XY = 2k,04

if hitran_version == 'h92'
  linker = '92';
elseif hitran_version == 'h96'
  linker = '96';
elseif hitran_version == 'h98'
  linker = '98';
elseif hitran_version == 'h2k'
  linker = '00';
elseif hitran_version == 'h04'
  linker = '04';
elseif hitran_version == 'h08'
  linker = '08';
  end

symlinker = ['!/bin/ln -s qtips' linker '.m qtips.m']; 
linker0 = '98';

ee = exist('qtips.m','file');
if (ee == 0 & (length(intersect(hitran_version,{'h92','h96','h98'})) == 1))
  symlinker = ['!/bin/ln -s qtips' linker0 '.m qtips.m']; 
  fprintf(1,'qtips.m not found .. making VERS %s symbolic link! \n',linker0);
  eval(symlinker);
elseif (ee ~= 0 & (length(intersect(hitran_version,{'h92','h96','h98'})) == 1))
  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  lser = ['!ls -lt qtips.m >& ' frand];
  eval(lser)
  fid = fopen(frand,'r');
  ugher = fscanf(fid,'%s');
  fclose(fid);
  rmer = ['!/bin/rm ' frand];
  eval(rmer);
  eee = strcmp('lrwxrwxrwx',ugher(1:10));
  if eee == 0
    error('oh oh qtips.m is NOT a symbolic link; do not want to delete it!!!')
  else
    %%check to see if the symbolic link is for the right version
    tata = findstr('->qtips',ugher);
    oldvers = ugher(tata+7:tata+8);
    fff = strcmp(oldvers,linker);
    %% if fff == 1, the old version == current wanted version ==> do nothing!
    if fff == 0
       %% symbolic link is to a different version than current wanted version
      fprintf(1,'updating VERS %s symbolic link for qtips.m ! \n',linker)
      rmer = ['!/bin/rm qtips.m']; eval(rmer);
      %fprintf(1,'symbolic link command : %s \n',symlinker);
      eval(symlinker);
      end
    end
elseif (ee ~= 0 & (length(intersect(hitran_version,{'h92','h96','h98'})) ~= 1))
  %% if we are running H2k or H04, just make sure you update to qtips04.m
  linker0 = '04';
  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  lser = ['!ls -lt qtips.m >& ' frand];
  eval(lser)
  fid = fopen(frand,'r');
  ugher = fscanf(fid,'%s');
  fclose(fid);
  rmer = ['!/bin/rm ' frand];
  eval(rmer);
  eee = strcmp('lrwxrwxrwx',ugher(1:10));
  if eee == 0
    error('oh oh qtips.m is NOT a symbolic link; do not want to delete it!!!')
  else
    %%check to see if the symbolic link is for the right version
    tata = findstr('->qtips',ugher);
    oldvers = ugher(tata+7:tata+8);
    fff = strcmp(oldvers,linker0);
    %% if fff == 1, the old version == current wanted version ==> do nothing!
    if fff == 0
       %% symbolic link is to a different version than current wanted version
      symlinker = ['!/bin/ln -s qtips' linker0 '.m qtips.m']; 
      fprintf(1,'just making VERS %s symbolic link for qtips.m ! \n',linker0)
      rmer = ['!/bin/rm qtips.m']; eval(rmer);
      eval(symlinker);
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
