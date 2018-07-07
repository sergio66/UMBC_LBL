function [kx, freq, temp] = contread(fname);

% function [kx, freq, temp] = contread(fname);
%
% Read a kCARTA continuum lookup table
%
% Input:
%    fname = [string] name of file to read
%
% Output:
%    kx    = [ntemp x nfreq] optical depth
%    freq  = [    1 x nfreq] frequency
%    temp  = [ntemp x     1] temperature
%
% comment: the continuum database files are in dir
%    /asl/data/kcarta/KCARTADATA/General/CKDieee_le
%

% Created: 28 May 2003 Howard Motteler
% Update: 12 Jun 2009, Scott Hannon - replace "mdeal" function with
%    reads of individual varibles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fname, 'r');
if fid == -1
  error(sprintf('can not open %s', fname));
end

%%header = f1,f2,df
% filemark=8+8+8;
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,[fstart,fend,df],'real*8');
% fwrite(fid,filemark,'integer*4');

%j = fread(fid, 1, 'integer*4');
%t = fread(fid, 3, 'real*8');
%j = fread(fid, 1, 'integer*4');
%[fstart, fend, df] = mdeal(t);
j      = fread(fid, 1, 'integer*4');
fstart = fread(fid, 1, 'real*8');
fend   = fread(fid, 1, 'real*8');
df     = fread(fid, 1, 'real*8');
j      = fread(fid, 1, 'integer*4');


%%header = loc,CKD vers,m,n
% filemark=4+4+4+4;
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,[loc,CKD,m,n],'integer*4');
% fwrite(fid,filemark,'integer*4');

%j = fread(fid, 1, 'integer*4');
%t = fread(fid, 4, 'integer*4');
%j = fread(fid, 1, 'integer*4');
%[loc, CKD, m, n] = mdeal(t);
j   = fread(fid, 1, 'integer*4');
loc = fread(fid, 1, 'integer*4');
CKD = fread(fid, 1, 'integer*4');
m   = fread(fid, 1, 'integer*4');
n   = fread(fid, 1, 'integer*4');
j   = fread(fid, 1, 'integer*4');


%save the temps
% temp=100:10:400;
% filemark = 8*length(temp);
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,temp,'real*8');
% fwrite(fid,filemark,'integer*4');

j = fread(fid, 1, 'integer*4');
ntemp = round(j/8); % exact integer
temp = fread(fid, ntemp, 'real*8');
j = fread(fid, 1, 'integer*4');

%%save the matrix 
%%%%loop over the rows (temperature)
% filemark = 8*length(fr);
% for i = 1:m
%    fwrite(fid,filemark,'integer*4');
%    fwrite(fid,ks(i,:)','real*8'); %'
%    fwrite(fid,filemark,'integer*4');
% end

kx = zeros(m,n);

for i = 1 : m
  j = fread(fid, 1, 'integer*4');
  kx(i,:) = fread(fid, [1,n], 'real*8');
  j = fread(fid, 1, 'integer*4');
end

fclose(fid);

%freq = fstart : df : fend;

freq = 1:n;
freq = fstart +(freq-1)*df;

%[min(freq) fstart max(freq) fend]
%[min(freq)/fstart max(freq)/fend]

%%% end of function %%%
