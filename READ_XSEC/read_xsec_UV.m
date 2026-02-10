function xs = read_xsec_UV(gid, gd, f1, f2)

% function xs = read_xsec_UV(gid, gd)
%
% read_xsec_UV() reads a HITRAN UV cross-section data file
% and returns a matlab structure containing the data; it
% assumes one gas per file.  
%
% input
%   gid  - HITRAN xsec gas file or gas ID
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
% read_xsec_UV() can deal with minor variations in file formats found 
% in the HITRAN database; in particular, it will skip over trailing
% zeros fields at the end of the last absorption record, and doesn't
% mind DOS format files.  It assumes that the first six header fields
% are 10-character records, and that header field five contains
% meaningful pressure info.  (The pressure can be zero, but no text
% should appear in this field.)
% 

% H. Motteler,  3 Apr 00
% S. Machado,   3 Jan 10
% Avogadro's number * 10^3
avogkm = 6.0221367e+26;

if gid == 7
  disp('this is gasID 7 ==> HITRAN units are in (cm2/molecule)^2')
  avogkm = avogkm^2;
  end

% set default directory for xsec data
%if nargin == 1
  gd = [hitranpath '/H2024/IR-XSect/Compressed-files/'];              idd = 2024;

  gd = '/asl/data/hitran/HITRAN04/UV/Xsec';
  gd = '/spinach/s6/sergio/RUN8_NIRDATABASE/UV/Xsec/';  
  gd = [hitranpath '/HITRAN04/UV/XSec/'];                             idd = 2004;                    
%end

% if given numeric gas ID, translate to gas filename
% see /home/sergio/SPECTRA/findxsec_plot_UV.m
if isnumeric(gid)
  switch gid
  % gas IDs from GENLN2
    case 3,   gstr = 'O3';
    case 7,   gstr = 'O2-O2';
    case 9,   gstr = 'SO2';
    case 10,  gstr = 'NO2';
    case 20,  gstr = 'H2CO'; 
    otherwise, error('unknown gas id')
    end
  end

if gid == 7
  if f2 <= 14000
    gstr = [gstr '_nir'];
  elseif f1 >= 14000
    gstr = [gstr '_vis'];  
    end
  end

if gid == 7
  fnamePOST='_UV04.xsc';
elseif gid == 3 | gid == 20
  fnamePOST='-UV04.xsc';
elseif gid == 9 & f1 >= 41691.05 
  fnamePOST='-UV00.xsc'; 
%elseif gid == 9 & f1 >= 23995.000 
elseif gid == 9
  fnamePOST='-UV08-UVA.xsc'; 
elseif gid == 10 
  fnamePOST='-UV00.xsc'; 
else
  gid
  error('huh ... no fnamePOST');
  end
gfname = [gstr fnamePOST];

% open the xsec gas data file
gfname = [gd, '/', gfname];
[fid, msg] = fopen(gfname, 'r');
if fid == -1
  fprintf(1,'whoops error opening file %s \n',gfname)
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
  gstr = xhead(1:10);
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
  if v1 < 8500 | v2 > 65000 | temp < 100 | 400 < temp 
    v1, v2, temp
    error('bad xsec header record');
  end
    
  % read absorption data
  [absc,j] = fscanf(fid, '%g', npts);
  if j ~= npts
    [j npts]
    error('short read of absorption data');
  end

  % get next header record
  % we may have to skip over some trailing junk, so use the 
  % initial gas id string to indicate a valid header record,
  % and keep reading until we see this or reach end of file
  xhead = fgetl(fid);
  while isstr(xhead) & (length(xhead) < 10 | ~strcmp(gstr, xhead(1:10)) )
    xhead = fgetl(fid);
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
  if gid ~= 7
    xs(nrec, band).absc = absc * avogkm;
  else
    disp('WARNING : Gas 7 : abscf in units read_xsec_UV')
    %% xs(nrec, band).absc = sqrt(max(absc,0)) * avogkm;
    xs(nrec, band).absc = absc * avogkm;
    end
  end
fclose(fid);

