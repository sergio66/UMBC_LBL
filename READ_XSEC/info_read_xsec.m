function [data,info] = info_read_xsec(fname);

% function [data,info] = info_read_xsec(fname);
%
% Read HITRAN 2008 IR cross-section data
%
% Input:
%    fname  -  [string] name of xsec file to read
%
% Output:
%    data   -  [m x n] absorption {cm^2/molecules}
%    info   -  [structure] supplemental info with fields:
%       chemical_symbol - [string] chemical symbol
%       common_name     - [string] common name
%       broadening_gas  - [string] broadening gas
%       filename        - [string] fname (without path)
%       reference_num   - [1 x 1] reference number
%       npts            - [1 x n] number of points, npts <= m
%       fmin            - [1 x n] min freq {cm^-1}
%       fmax            - [1 x n] max freq {cm^-1}
%       pressure        - [1 x n] pressure {torr}
%       temperature     - [1 x n] temperature {Kelvin}
%       FTS_resolution  - [1 x n] resolution of FTS measurement {cm^-1}
%       max_xsec        - [1 x n] maximum value of data
%

% Created: 02 July 2010, Scott Hannon
% Update: 14 Jul 2010, S.Hannon - bug fix for ncount > npoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxm = 60000;
maxn = 55;


% Check input
d = dir(fname);
if (length(d) ~= 1)
   error(['Unable to open fname: ' fname])
end
info.filename = d.name;


% Declare output
data = zeros(maxm, maxn);
info.npts = zeros(1,maxn);
info.fmin = zeros(1,maxn);
info.fmax = zeros(1,maxn);
info.pressure = zeros(1,maxn);
info.temperature = zeros(1,maxn);
info.FTS_resolution = zeros(1,maxn);
info.max_xsec = zeros(1,maxn);


% Open file
fid = fopen(fname,'r');

lheader = 1; % header
ncount = 0;
% Loop over rows in file
while 1
   tline = fgetl(fid); % read current line
   if ~ischar(tline), break, end  % abort if not a string
   if (lheader == 1)
      ncount = ncount + 1; % new table
      % read header
tline
      info.chemical_symbol = tline(1:20);
      info.fmin(ncount) = str2num( tline(21:30) );
      info.fmax(ncount) = str2num( tline(31:40) );
      info.npts(ncount) = str2num( tline(41:47) );
      info.temperature(ncount) = str2num( tline(48:54) );
      info.pressure(ncount) = str2num( tline(55:60) );
      info.max_xsec(ncount) = str2num( tline(61:70) );
      info.FTS_resolution(ncount) = str2num( tline(71:75) );
      info.common_name = tline(76:90);
      % Note: tline(91-94) not used
      info.broadening_gas = tline(95:97);
      info.reference_num = str2num( tline(98:100) );
      %
      lheader = 0;
      mcount = 0;
      %
   else
      % read data
      [junk, mnew] = sscanf(tline,'%e');
      ind = mcount + (1:mnew);
      data(ind,ncount) = junk;
      mcount = round(mcount + mnew); % exact integer
      if (mcount >= info.npts(ncount))
         % Done with current table; next line is header
         lheader = 1;
      end
   end

disp(['ncount,mcount = ' int2str(ncount) ', ' int2str(mcount)])

end
fclose(fid);


% Shrink pre-declared output to size actually required
amaxm = max(info.npts);
amaxn = ncount;
indm = 1:amaxm;
indn = 1:amaxn;
data = data(indm,indn);
info.npts = info.npts(indn);
info.fmin = info.fmin(indn);
info.fmax = info.fmax(indn);
info.pressure = info.pressure(indn);
info.temperature = info.temperature(indn);
info.FTS_resolution = info.FTS_resolution(indn);
info.max_xsec = info.max_xsec(indn);

%%% end of function %%%
