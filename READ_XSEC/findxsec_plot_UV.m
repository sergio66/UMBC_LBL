function [iYes,xsecfilename] = findxsec_plot_UV(wv1,wv2,gid);

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
tv = ttorad(wv,t0);

gasID = gid;

fnamePRE = '/asl/data/hitran/HITRAN04/UV/Xsec/';
fnamePRE = '/spinach/s6/sergio/RUN8_NIRDATABASE/UV/Xsec/';
fnamePRE = [hitranpath '/HITRAN04/UV/Xsec/'];

%% from ~sergio/KCARTA/DOC/gasids
%% also see /home/sergio/abscmp/READ_XSEC/gid2mol.m
%gname{1} = 'Br0';   % 
%gname{2} = 'H2CO';  % 20
%gname{3} = 'NO3';   % 
%gname{4} = 'O2-O2'; % 
%gname{5} = 'O3';    % 3
%gname{6} = 'OClO';  % 

switch gid
% gas IDs from GENLN2
  case 3,  gstr = 'O3';
  case 7,  gstr = 'O2-O2';
  case 9,  gstr = 'SO2';
  case 10,  gstr = 'NO2';
  case 20, gstr = 'H2CO'; 
  otherwise, error('unknown gas id')
end

if gid == 7
  if wv2 <= 14000
    gstr = [gstr '_nir'];
  elseif wv1 >= 14000
    gstr = [gstr '_vis'];  
    end
  end
    
fnameIN = gstr;

if gid == 3
  %% this gives cross sections in units of cm2/molecule
  hitlin_fname = '/asl/data/hitran/HITRAN04/UV/Xsec/Alt/O3-UV04.alt';
  %% /asl/data/hitran/HITRAN04/UV/Xsec/Alt/O3-UV04.alt says cross section
  %% at 245.107 nm is sigmaHIT (cm2/molecule) = 1e-20 * 9.9079E2 at 273 K 

  %% supplied by Li Zhu 2008, gives cross sections in units of /(atm cm)
  hitlin_fname = '/home/sergio/SPECTRA/VISIBLE_OD/O3.csv';
  %% recall pV = nRT ==> n/V = 2.4486e+19 molecules/cm3 at 296 K
  %%        pV = nRT ==> n/V = 2.6549e+19 molecules/cm3 at 273 K
  %% so this file says cross section sigmaLI = 267.063 (1/atm cm) at 245.1 nm
  %% sigmaLI(1/atm cm) = sigmaHIT(cm2/molecule) * (molecules/cm3)/atm
  %%                   = 1e-20 * 9.9079E2 * 2.6549e+19

  liO3 = load(hitlin_fname);
  liO3(:,1) = liO3(:,1)/10;     %% change from Angstroms to nm
  liO3(:,1) = liO3(:,1)/1000;   %% change from nm to um
  liO3(:,1) = 10000./liO3(:,1); %% change from um to wavenumber 
  v1 = min(liO3(:,1));
  v2 = max(liO3(:,1));

elseif gid ~= 3
  if gid == 7
    fnamePOST='_UV04.xsc';
  elseif gid == 3 | gid == 20
    fnamePOST='-UV04.xsc';
  elseif gid == 9 & wv1 >= 41691.05
    fnamePOST='-UV00.xsc';
%  elseif gid == 9 & wv1 >= 23995.000
  elseif gid == 9
    fnamePOST='-UV08-UVA.xsc';
  elseif gid == 10
    fnamePOST='-UV00.xsc';
  else
    gid
    error('huh ... no fnamePOST');
    end
  hitlin_fname=[fnamePRE fnameIN fnamePOST];

  ee = exist(hitlin_fname);
  if ee == 0
    fprintf(1,'looking for %s %3i \n',hitlin_fname,ee);
    error('could not find file');
    end

  start = wv1;
  stop  = wv2;

  iYes = -1;

  clf
  %header = ['!head -1 ' hitlin_fname]; eval(header);

  %% reader is copied from abscmp/READ_XSEC/read_xsec.m
  gf = hitlin_fname;
  [fid, msg] = fopen(gf, 'r');
  if fid == -1
    error(msg);
  end

  % read first header
  xhead = fgetl(fid); 

  k = length(xhead);

  % get xhead fields  --- these are slightly different from the IR xsec header
  gxstr= xhead(1:10);
  v1   = str2num(xhead(21:30));
  v2   = str2num(xhead(31:40));
  npts = str2num(xhead(41:48));
  temp = str2num(xhead(49:54));
  pres = str2num(xhead(55:61));
  if length(pres) == 0
    pres = 0;   %% case for gid 7 has "varies"
    end
  maxi = str2num(xhead(62:70));
  junk = xhead(71:k);

  % sanity check
  if v1 < 8500 | v2 > 52000 | temp < 100 | 400 < temp 
    v1, v2, temp
    error('bad xsec header record');
    end
  fclose(fid);
  end

iYes = -1;
if wv1 <= v2 & wv2 >= v1
  iYes = 1;
  end
%[v1 v2 iYes]

xsecfilename = gstr;
