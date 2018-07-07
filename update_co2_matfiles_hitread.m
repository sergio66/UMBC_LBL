iLinkBands = -1;
iCO2       = -1;
if gasID == 2
  iCO2 = 1;
end

if iCO2 > 0  
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
    error('sorry .. could not figure out which dir to link CO2_FILES');
  end

  randstr = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  frand   = ['ugh' randstr];
  frand   = mktemp('ugh');

  CO2matfilesdir0 = pwd; 
  CO2matfilesdir = pwd; 
  CO2matfilesdir = [CO2matfilesdir '/CO2_MATFILES/'];

  cder = ['cd ' CO2matfilesdir]; eval(cder);
  lser = ['!ls -lt hit_info2350.mat >& ' frand];
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
  cder = ['cd ' CO2matfilesdir0]; eval(cder);
end

if (iLinkQTIPS == 1 | iLinkMASS == 1 | iLinkBands == 1) & (iCO2 > 0)
  disp('need to link the CO2 band files in /home/sergio/SPECTRA/CO2_MATFILES to correct HXY version')
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
    error('sorry .. could not figure out which dir to link CO2_FILES');
  end

  % may as well update the symlinks for CO2_MATFILES
  disp('----> updating CO2_MATFILES symbolic links')
  CO2matfilesdir0 = pwd; 
  CO2matfilesdir = pwd; 
  CO2matfilesdir = [CO2matfilesdir '/CO2_MATFILES/'];
  tempFname = [CO2matfilesdir 'link_hitco2.txt'];
  if ~exist(tempFname)
    fprintf(1,'% DNE \n',tempFname)
    error('ohoh')
  end    
  glist = textread(tempFname,'%s','delimiter','\n','whitespace','');
  cder = ['cd ' CO2matfilesdir]; eval(cder);
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
  cder = ['cd ' CO2matfilesdir0]; eval(cder);
end

cder = ['cd ' current_dir]; eval(cder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
this is CO2_MATFILES/link_hitco2.txt

hit1932.mat
hit2080.mat
hit2093.mat
hit2129.mat
hit2310.mat
hit2311.mat
hit2320.mat
hit2321.mat
hit2322.mat
hit2350.mat
hit2351.mat
hit2352.mat
hit2353.mat
hit2354.mat
hit618.mat
hit648.mat
hit662.mat
hit667.mat
hit668.mat
hit720.mat
hit740.mat
hit791.mat
hit_info1932.mat
hit_info2080.mat
hit_info2093.mat
hit_info2129.mat
hit_info2310.mat
hit_info2311.mat
hit_info2320.mat
hit_info2321.mat
hit_info2322.mat
hit_info2350.mat
hit_info2351.mat
hit_info2352.mat
hit_info2353.mat
hit_info2354.mat
hit_info618.mat
hit_info648.mat
hit_info662.mat
hit_info667.mat
hit_info668.mat
hit_info720.mat
hit_info740.mat
hit_info791.mat

%}
