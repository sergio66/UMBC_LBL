iLinkBands = -1;
iCO       = -1;
if gasID == 5
  iCO = 1;
end

if iCO > 0  
  if strcmp('h92',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h96',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h98',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h2k',hitran_version) 
    HITdir = 'H2K/';
  elseif strcmp('h04',hitran_version) 
    HITdir = 'H04/';
  elseif strcmp('h08',hitran_version) 
    HITdir = 'H08/';
  elseif strcmp('h12',hitran_version) 
    HITdir = 'H12/';
  elseif strcmp('h16',hitran_version) 
    HITdir = 'H16/';
  else
    error('sorry .. could not figure out which dir to link CO_FILES');
  end

  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  frand   = mktemp('ugh');

  COmatfilesdir0 = pwd; 
  COmatfilesdir = pwd; 
  COmatfilesdir = [COmatfilesdir '/CO_MATFILES/'];

  cder = ['cd ' COmatfilesdir]; eval(cder);
  lser = ['!ls -lt hit_info2150.mat >& ' frand];
  eval(lser)
  fid = fopen(frand,'r');
  ugher = fscanf(fid,'%s');
  fclose(fid);
  rmer = ['!/bin/rm ' frand];
  eval(rmer);

  tata = findstr('->',ugher);
  oldvers = ugher(tata+2:tata+5);
%  fprintf(1,'debugging hitread.m .... bands linker = %s oldvers = %s need to be linked to %s \n',ugher,oldvers,HITdir)
  fff = strcmp(oldvers,HITdir);
  if fff > 0
    %% oldvers and neededvers agree
    iLinkBands = -1;
  else
    %% oldvers and neededvers disagree  
    iLinkBands = +1;
  end
  cder = ['cd ' COmatfilesdir0]; eval(cder);
end

if (iLinkQTIPS == 1 | iLinkMASS == 1 | iLinkBands == 1) & (iCO > 0)
  disp('need to link the CO band files in /home/sergio/SPECTRA/CO_MATFILES to correct HXY version')
  if strcmp('h92',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h96',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h98',hitran_version) 
    HITdir = 'H98/';
  elseif strcmp('h2k',hitran_version) 
    HITdir = 'H2K/';
  elseif strcmp('h04',hitran_version) 
    HITdir = 'H04/';
  elseif strcmp('h08',hitran_version) 
    HITdir = 'H08/';
  elseif strcmp('h12',hitran_version) 
    HITdir = 'H12/';
  elseif strcmp('h16',hitran_version) 
    HITdir = 'H16/';
  else
    error('sorry .. could not figure out which dir to link CO_FILES');
  end

  % may as well update the symlinks for CO_MATFILES
  disp('----> updating CO_MATFILES symbolic links')
  COmatfilesdir0 = pwd; 
  COmatfilesdir = pwd; 
  COmatfilesdir = [COmatfilesdir '/CO_MATFILES/'];
  tempFname = [COmatfilesdir 'link_hitco.txt'];
  if ~exist(tempFname)
    fprintf(1,'% DNE \n',tempFname)
    error('ohoh')
  end    
  glist = textread(tempFname,'%s','delimiter','\n','whitespace','');
  cder = ['cd ' COmatfilesdir]; eval(cder);
  for ii = 1 : length(glist)
    if exist(glist{ii}) & exist([HITdir glist{ii}])
      rmer = ['!/bin/rm ' glist{ii}];
      lner = ['!ln -s   ' HITdir glist{ii}  ' ' glist{ii}];
      eval(rmer);
      eval(lner);
    elseif ~exist(glist{ii}) & exist([HITdir glist{ii}])
      lner = ['!ln -s   ' HITdir glist{ii}  ' ' glist{ii}];
      eval(lner);
    elseif exist(glist{ii}) & ~exist([HITdir glist{ii}])
      fprintf(1,'%s EXISTS %s DNE \n',glist{ii},[HITdir glist{ii}])    
    elseif ~exist(glist{ii}) & ~exist([HITdir glist{ii}])
      fprintf(1,'%s or %s DNE \n',glist{ii},[HITdir glist{ii}])
    end
  end
  cder = ['cd ' COmatfilesdir0]; eval(cder);
end

cder = ['cd ' current_dir]; eval(cder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
this is CO_MATFILES/link_hitco.txt

hit2150.mat
hit_info2150.mat

%}
