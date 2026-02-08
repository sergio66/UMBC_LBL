function xs = read_xsec(gf, gd, HITRAN)

iPrint = -1;

% function xs = read_xsec(gf, gd, HITRAN)
%
% read_xsec() reads a HITRAN IR cross-section data file
% and returns a matlab structure containing the data; it
% assumes one gas per file.  
%
% input
%   gf  - HITRAN xsec gas file or gas ID
%   gd  - optional directory for gas files 
%   HITRAN (optional) = 2008 or 20012 or 2016 (default)
%
% output
%   xs  - an nrec by nband structure of xsec records
% 
% xs is an nrec by nband array of structures.  The structure fields
% are the xsec record header values and the absorption data for that
% record.  The fields returned are:
% 
%   gstr  - gas string identifier
%   v1    - band starting wavenumber
%   v2    - band ending wavenumber
%   npts  - number of frequency points
%   temp  - temperature
%   pres  - pressure
%   absc  - a npts x 1 array of absorption data, kmoles/cm^2
%
% read_xsec() can deal with minor variations in file formats found 
% in the HITRAN database; in particular, it will skip over trailing
% zeros fields at the end of the last absorption record, and doesn't
% mind DOS format files.  It assumes that the first six header fields
% are 10-character records, and that header field five contains
% meaningful pressure info.  (The pressure can be zero, but no text
% should appear in this field.)
% 

% H. Motteler,  3 Apr 00

% Avogadro's number * 10^3
avogkm = 6.0221367e+26;

if nargin <= 2
  HITRAN = 2012;
  HITRAN = 2016;  
  HITRAN = 2020;
  HITRAN = 2024;    
end

% set default directory for xsec data
%if nargin == 1
  if HITRAN == 1998
    gd = '/asl/data/hitran/xsec98.ok';  idd = 1998;
  elseif HITRAN == 2008
    gd = '/asl/data/hitran/HITRAN2008/IR-XSect/Compressed-files/Junk/'; idd=2008;
    gd = '/asl/data/hitran/HITRAN08_SERGIO/Xsec/'; idd=2008;
    gd = '/asl/data/hitran/H2008/IR-XSect/Uncompressed-files/'; idd=2008;
  elseif HITRAN == 2012 
    gd = '/asl/data/hitran/H2012/IR-XSect/Uncompressed-files/'; idd=2012;
  elseif HITRAN == 2015 
    gd = '/asl/data/geisa/G2015/2015.IR-XSect/Uncompressed-files/'; idd=2015;
  elseif HITRAN == 2016 
    gd = '/asl/data/hitran/H2016/IR-XSect/Uncompressed-files/'; idd=2016;
  elseif HITRAN == 2020 
    gd = '/asl/data/hitran/H2020/IR-XSect/Uncompressed-files/'; idd=2020;
    gd = '/umbc/xfs3/strow/asl/rta/hitran/H2020/IR-XSect/Uncompressed-files/'; idd=2020;        
  elseif HITRAN == 2024 
    gd = '/asl/data/hitran/H2020/IR-XSect/Uncompressed-files/'; idd=2024;
    gd = '/umbc/xfs3/strow/asl/rta/hitran/H2024/IR-XSect/Uncompressed-files/'; idd=2024;    
  else
    error('need HITRAN == 1998 2008 2012 (G)2015 2016 2020')
  end
%end

% if given numeric gas ID, translate to gas filename
if isnumeric(gf)
  fprintf(1,'  read_xsec gid = %3i HITRAN = %4i \n',gf,HITRAN)
  gf = sprintf('%s.xsc', gid2mol(gf,HITRAN));
end

% open the xsec gas data file
gf = [gd, '/', gf];
fprintf(1,'looking to open and read %s \n',gf)
if ~exist(gf)
  fprintf(1,'OOPS %s does not exist !!! Exiting read_xsec.m!!! \n',gf);
  xs = [];
  return
end

[fid, msg] = fopen(gf, 'r');
if fid == -1
  error(msg);
end

nrec = 0;
band = 0;
v1prev = 0;
v2prev = 0;

% read first header
xhead = fgetl(fid); 

while isstr(xhead)

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
    %% see hitran-xsec.pdf or http://hitran.org/docs/cross-sections-definitions/
    %% example
    %%                CH3CN 3881.8386 4573.9200  11484  276.1 760.0 1.162E-20 0.11   Acetonitrile     N2 32
    %% 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    %%          1         2         3         4         5         6         7         8         9        10
    %%
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
  
  % sanity check
  if v1 < 10 | 50000 <  v2 | temp < 100 | 400 < temp 
    v1, v2, temp
    error('bad xsec header record');
  end
    
  % read absorption data
  [absc,j] = fscanf(fid, '%g', npts);
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
  elseif idd > 1998
    while isstr(xhead) & (length(xhead) < 20 | ~strcmp(gstr, xhead(1:20)) )
      xhead = fgetl(fid);
    end
  end

  % see if we're in a new band
  %% if v2prev < v1   ORIG by Howard Motteler    
  if round(v2prev) < round(v1) | round(v1prev) < round(v1)
    nrec = 0;
    band = band + 1;
  end

  v1prev = v1;
  v2prev = v2;
  
  % increment record count for this band
  nrec = nrec + 1;

  if iPrint > 0
    fprintf(1,'band %3i nrec = %3i v1,v2 = %8.3f %8.3f cm-1, npts=%8i temp=%8.3f K pres=%8.3f mb \n',band,nrec,v1,v2,npts,temp,pres)
  end
  
  % save this record in returned structure
  xs(nrec, band).gstr = gstr;
  xs(nrec, band).v1   = v1;
  xs(nrec, band).v2   = v2;
  xs(nrec, band).npts = npts;
  xs(nrec, band).temp = temp;
  xs(nrec, band).pres = pres;
  xs(nrec, band).maxi = maxi;
  xs(nrec, band).absc = absc * avogkm;

end
fclose(fid);
