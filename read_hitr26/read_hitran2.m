
% NAME
%   read_hitran2 - read HITRAN data, matlab version
%
% SYNOPSIS
%   s = read_hitran2(v1, v2, str, gid, hsrc);
%
% INPUTS
%   v1   - wavenumber lower bound
%   v2   - wavenumber upper bound
%   str  - line strength lower bound
%   gid  - HITRAN gas id 
%   hsrc - HITRAN data file or directory
%
% OUTPUTS
%   s.igas    - molecule ID number
%   s.iso     - isotope number
%   s.wnum    - line wavenumber
%   s.stren   - line strength
% 
%   s.tprob   - transition probability
%   s.abroad  - air-broad half width
%   s.sbroad  - self-broad half width
%   s.els     - lower-state energy
%   s.abcoef  - coeff of temp dep of ABROAD
% 
%   s.tsp     - transition shift due to pressure 
%   s.iusgq   - upper state global quanta index
%   s.ilsgq   - lower state global quanta index
%   s.uslq    - upper state local quanta
%   s.bslq    - lower state local quanta
%   s.ai      - accuracy indices
%   s.ref     - indices for lookup of references
%   s.flag    - flag
%   s.swus    - statistical weight, upper state
%   s.swls    - statistical weight, lower state
% 
% NOTES
%   The optional input parameter hsrc can be either a HITRAN data file,
%   or a directory of HITRAN files organized gas types.  For the latter
%   case, the individual files are assumed to have have names g<n>.dat,
%   for gas <n>.  It is possible (but very slow) to specify the whole,
%   unsplit HITRAN database file as hsrc, and select gasses by gid; the
%   usual way to use the procedure is to specify a directory for hsrc,
%   and select gasses by gasid.  
% 
%   The default name for hsrc is "hitran.dat"; this is convenient for
%   setting the gas directory with a symlink.
%
% AUTHOR
%   H. Motteler, 15 Dec 06
%

function s = read_hitran2(v1, v2, str, gid, hsrc);

if nargin == 4
 hsrc = 'hitran.dat';
end

if exist(hsrc) == 7
  hsrc = [hsrc, '/g', num2str(gid), '.dat'];
end

[fid,msg] = fopen(hsrc, 'r');
if fid == -1
  error(msg)
end

[A, k] = fscanf(fid, '%161c', [161, inf]);

A = A';

igas  = str2num(A(:, 1:2));
wnum  = str2num(A(:, 4:15));
stren = str2num(A(:, 16:25));

q = (igas == gid) & (stren >= str) & (v1 <= wnum) & (wnum <= v2);

A = A(q,:);

s.igas    =   str2num(A(:, 1:2));
s.iso     =   str2num(A(:, 3:3));
s.wnum    =   str2num(A(:, 4:15));
s.stren   =   str2num(A(:, 16:25));

s.tprob   =   str2num(A(:, 26:35));
s.abroad  =   str2num(A(:, 36:40));
s.sbroad  =   str2num(A(:, 41:45));
s.els     =   str2num(A(:, 46:55));
s.abcoef  =   str2num(A(:, 56:59));

s.tsp     =   str2num(A(:, 60:67));
s.iusgq   =   A(:, 68:82);
s.ilsgq   =   A(:, 83:97);
s.uslq    =   A(:, 98:112);
s.bslq    =   A(:, 113:127);
s.ai      =   A(:, 128:133);
s.ref     =   A(:, 134:145);
s.flag    =   A(:, 146:146); 
s.swus    =   str2num(A(:, 147:153));
s.swls    =   str2num(A(:, 154:160));

