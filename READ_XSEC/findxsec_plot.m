function [iYes,xsecfilename,outhead] = findxsec_plot(wv1,wv2,gid,HITRANyear);

xhead = struct;

if nargin < 4
  HITRANyear = 2012;
  HITRANyear = 2016;
  HITRANyear = 2020;
end

%addpath /home/sergio/HITRAN2UMBCLBL/READ_XSEC
addpath /home/sergio/SPECTRA/READ_XSEC
%addpath /home/sergio/SPECTRA/READ_XSEC_GEISA

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
tv = ttorad(wv,t0);

%% see /home/sergio/abscmp/READ_XSEC/gid2mol.m
%% gname{id} = gid2mol(gid,HITRANyear);

id = gid-51+1;
if id >= 1 & id <= 31 & mod(HITRANyear,4) == 0
  disp('looking for HITRAN xcec')
  fnameIN = gid2mol(gid,HITRANyear);
elseif id >= 1 & id <= 31 & mod(HITRAN,4) ~= 4
  disp('looking for GEISA xcec')
  fnameIN = gid2mol(gid,HITRANyear);  
else
  error('invalid xsec gasid : must be between 51-81');
end

%% see /home/sergio/HITRAN2UMBCLBL/READ_XSEC/read_xsec.m

if HITRANyear == 1998
  fnamePRE ='/asl/data/hitran/xsec98.ok/'; idd = 1998;
  fnamePRE ='/umbc/xfs3/strow/asl/rta/hitran/xsec98.ok/'; idd = 1998;
  fnamePRE = [hitranpath '/xsec98.ok/']; idd = 1998;
elseif HITRANyear == 2008
  fnamePRE = '/asl/data/hitran/HITRAN2008/IR-XSect/Compressed-files/Junk/'; idd = 2008;
  fnamePRE = '/asl/data/hitran/HITRAN08_SERGIO/Xsec/';                      idd = 2008;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/HITRAN08_SERGIO/Xsec/';       idd = 2008;
  fnamePRE = [hitranpath '/HITRAN08_SERGIO/Xsec/'];                         idd = 2008;  
elseif HITRANyear == 2012
  fnamePRE = '/asl/data/hitran/H2012/IR-XSect/Uncompressed-files/';                idd = 2012;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2012/IR-XSect/Uncompressed-files/'; idd = 2012;
  fnamePRE = [hitranpath '/H2012/IR-XSect/Uncompressed-files/'];                   idd = 2012;    
elseif HITRANyear == 2016
  fnamePRE = '/asl/data/hitran/H2016/IR-XSect/Uncompressed-files/';                idd = 2016;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2016/IR-XSect/Uncompressed-files/'; idd = 2016;
  fnamePRE = [hitranpath '/H2016/IR-XSect/Uncompressed-files/'];                   idd = 2016;      
elseif HITRANyear == 2015
  fnamePRE = '/asl/data/geisa/G2015/2015.IR-XSect/Uncompressed-files/';            idd = 2015;   %% ?? idd = 2016?
  fnamePRE = '/asl/data/geisa/G2015/2015.IR-XSect/Uncompressed-files/';            idd = 2015;   %% ?? idd = 2016?
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/G2015/IR-XSect/Uncompressed-files/'; idd = 2015;   %% ?? idd = 2016?
  fnamePRE = [hitranpath '/G2015/IR-XSect/Uncompressed-files/'];                   idd = 2015;   %% ?? idd = 2016?
elseif HITRANyear == 2020
  fnamePRE = '/asl/data/hitran/H2020/IR-XSect/Uncompressed-files/';                idd = 2020;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2020/IR-XSect/Uncompressed-files/'; idd = 2020;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2020/IR-XSect/Uncompressed-files/'; idd = 2020;
  fnamePRE = [hitranpath '/H2020/IR-XSect/Uncompressed-files/'];                   idd = 2020;        
elseif HITRANyear == 2024
  fnamePRE = '/asl/data/hitran/H2024/IR-XSect/Uncompressed-files/';                idd = 2024;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2024/IR-XSect/Uncompressed-files/'; idd = 2024;
  fnamePRE = '/umbc/xfs3/strow/asl/rta/hitran/H2024/IR-XSect/Uncompressed-files/'; idd = 2024;
  fnamePRE = [hitranpath '/H2024/IR-XSect/Uncompressed-files/'];                   idd = 2024;        
end

fnamePOST='.xsc';
hitlin_fname=[fnamePRE fnameIN fnamePOST];

ee = exist(hitlin_fname);
if ee == 0
  fprintf(1,'looking for %s %3i \n',hitlin_fname,ee);
  error('could not find file');
else
  fprintf(1,'gid %3i using %s \n',gid,hitlin_fname);
end

start = wv1;
stop  = wv2;

clf
%header = ['!head -1 ' hitlin_fname]; eval(header);

%% reader is copied from abscmp/READ_XSEC/read_xsec.m
gf = hitlin_fname;
%fprintf(1,'in SPECTRA/findxsec_plot fname = %s \n',gf);
[fid, msg] = fopen(gf, 'r');
if fid == -1
  error(msg);
end

nrec = 0; 
band = 0; 
v2prev = 0; 

% read first header
xhead = fgetl(fid); 

iCnt = 0;
iYes = -1;
while isstr(xhead) 

  iCnt = iCnt + 1;
  k = length(xhead);

  % get xhead fields
  if idd <= 1998
    gstr = xhead(1:10);
    v1   = str2num(xhead(11:20));
    v2   = str2num(xhead(21:30));
    npts = str2num(xhead(31:40));
    temp = str2num(xhead(41:50));
    pres = str2num(xhead(51:60));
    maxi = str2num(xhead(61:70));
    junk = xhead(71:k);
  else
    gstr = xhead(1:20);
    v1   = str2num(xhead(21:30));
    v2   = str2num(xhead(31:40));
    npts = str2num(xhead(41:47));
    temp = str2num(xhead(48:54));
    pres = str2num(xhead(55:60));
    maxi = str2num(xhead(61:70));
    fts  = str2num(xhead(71:75));
    junk = xhead(76:k);
  end

  outhead.v1(iCnt) = v1;
  outhead.v2(iCnt) = v2;
  outhead.npts(iCnt) = npts;
  outhead.temp(iCnt) = temp;
  outhead.pres(iCnt) = pres;  

  % sanity check
  if v1 < 10 | 50000 <  v2 | temp < 100 | 400 < temp 
    v1, v2, temp
    error('bad xsec header record');
  end

  if wv1 <= v2 & wv2 >= v1
    iYes = 1;
  end
  %[v1 v2 iYes]

  % read absorption data 
  [absc,j] = fscanf(fid, '%g', npts);
  dv = (v2 - v1)/(npts-1);
  ind = 1:npts;
  ind = ind - 1;
  v = v1 + (ind)*dv;
  outhead.v{iCnt} = v;
  outhead.absc{iCnt} = absc;
  
  if j ~= npts 
    error('short read of absorption data'); 
  end 

  % get next header record 
  % we may have to skip over some trailing junk, so use the  
  % initial gas id string to indicate a valid header record, 
  % and keep reading until we see this or reach end of file 
  xhead = fgetl(fid); 
  if idd <= 1998
    while isstr(xhead) & (length(xhead) < 10 | ~strcmp(gstr, xhead(1:10)) ) 
      xhead = fgetl(fid); 
    end 
  elseif idd >= 1998
    while isstr(xhead) & (length(xhead) < 20 | ~strcmp(gstr, xhead(1:20)) ) 
      xhead = fgetl(fid); 
    end 
  end
 
  % see if we're in a new band 
  if v2prev < v1 
    nrec = 0; 
    band = band + 1; 
  end 

  v2prev = v2; 
   
  % increment record count for this band 
  nrec = nrec + 1; 

end

fclose(fid);
%xsecfilename = gname{id};
xsecfilename = fnameIN;

fprintf(1,'  have read %4i bands \n',nrec);
