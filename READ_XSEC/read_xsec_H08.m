function xs = read_xsec(gf, gd)

% function xs = read_xsec(gf, gd)
%
% read_xsec() reads a HITRAN IR cross-section data file
% and returns a matlab structure containing the data; it
% assumes one gas per file.  
%
% input
%   gf  - HITRAN xsec gas file or gas ID
%   gd  - optional directory for gas files 
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

% set default directory for xsec data
if nargin == 1
  gd = '/asl/data/hitran/xsec98.ok';  idd = 1998;
  gd = '/asl/data/hitran/HITRAN2008/IR-XSect/Compressed-files/Junk/'; idd=2008;
  gd = '/asl/data/hitran/HITRAN08_SERGIO/Xsec/'; idd=2008;
end

% if given numeric gas ID, translate to gas filename
if isnumeric(gf)
  gf = sprintf('%s.xsc', gid2mol(gf));
end

% open the xsec gas data file
gf = [gd, '/', gf];
fprintf(1,'looking to open and read %s \n',gf)
[fid, msg] = fopen(gf, 'r');
if fid == -1
  error(msg);
end

nrec = 0;
band = 0;
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
  if v2prev < v1
    nrec = 0;
    band = band + 1;
  end
  v2prev = v2;
  
  % increment record count for this band
  nrec = nrec + 1;

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
