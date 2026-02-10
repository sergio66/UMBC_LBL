function [iYes,xsecfilename] = findxsec_plot_OLD(wv1,wv2,gid);

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
tv = ttorad(wv,t0);

gasID = gid;

fnamePRE = '/asl/data/hitran/xsec98.ok/';
  fnamePRE ='/umbc/xfs3/strow/asl/rta/hitran/xsec98.ok/'; idd = 1998;
  fnamePRE = [hitranpath '/xsec98.ok/']; idd = 1998;

%% from ~sergio/KCARTA/DOC/gasids
%% also see /home/sergio/abscmp/READ_XSEC/gid2mol.m
gname{1} = 'CCl3F';   % 51  CFCl3   (F11) 
gname{2} = 'CCl2F2';  % 52          (F12)
gname{3} = 'CClF3';   % 53    (F13)
gname{4} = 'CF4';     % 54    (F14)
gname{5} = 'CHCl2F';  % 55    (F21)
gname{6} = 'CHClF2';  % 56    (F22)
gname{7} = 'C2Cl3F3'; % 57   (F113)
gname{8} = 'C2Cl2F4'; % 58 (F114)
gname{9} = 'C2ClF5';  % 59  (F115)
gname{10} = 'CCl4';   % 60
gname{11} = 'ClONO2'; % 61
gname{12} = 'N2O5';   % 62
gname{13} = 'HNO4';   % 63

id = gid-51+1;

%% also see /home/sergio/abscmp/READ_XSEC/gid2mol.m
%% gname{id} = gid2mol(gid);

if id >= 1 & id <= 13
  fnameIN = gname{id};
else
  error('invalid xsec gasid : must be between 51-63');
  end

fnamePOST='.xsc';
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

  % get xhead fields
  gstr = xhead(1:10);
  v1   = str2num(xhead(11:20));
  v2   = str2num(xhead(21:30));
  npts = str2num(xhead(31:40));
  temp = str2num(xhead(41:50));
  pres = str2num(xhead(51:60));
  maxi = str2num(xhead(61:70));
  junk = xhead(71:k);

  % sanity check
  if v1 < 10 | 10000 <  v2 | temp < 100 | 400 < temp 
    v1, v2, temp
    error('bad xsec header record');
  end

iYes = -1;
if wv1 <= v2 & wv2 >= v1
  iYes = 1;
  end
%[v1 v2 iYes]

fclose(fid);
xsecfilename = gname{id};
