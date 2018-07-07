function [kx, freq, temp] = cont_le2be(CKD,self_for);

%%% translates le to be files for CKD

if self_for == +1
  extension = [num2str(CKD) '.bin'];
  fname = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
  fname2= ['/asl/data/kcarta/KCARTADATA/General/CKDieee_be/CKDSelf' extension];
else
  extension = [num2str(CKD) '.bin'];
  fname = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor' extension];
  fname2= ['/asl/data/kcarta/KCARTADATA/General/CKDieee_be/CKDFor' extension];
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('reading in the file .....')
fid = fopen(fname, 'r');
if fid == -1
  error(sprintf('can not open %s', fname));
end

%%header = f1,f2,df
% filemark=8+8+8;
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,[fstart,fend,df],'real*8');
% fwrite(fid,filemark,'integer*4');

j = fread(fid, 1, 'integer*4');
t = fread(fid, 3, 'real*8');
j = fread(fid, 1, 'integer*4');
[fstart, fend, df] = mdeal(t);

%%header = loc,CKD vers,m,n
% filemark=4+4+4+4;
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,[loc,CKD,m,n],'integer*4');
% fwrite(fid,filemark,'integer*4');

j = fread(fid, 1, 'integer*4');
t = fread(fid, 4, 'integer*4');
j = fread(fid, 1, 'integer*4');
[loc, CKD, m, n] = mdeal(t);

%save the temps
% temp=100:10:400;
% filemark = 8*length(temp);
% fwrite(fid,filemark,'integer*4');
% fwrite(fid,temp','real*8');
% fwrite(fid,filemark,'integer*4');

j = fread(fid, 1, 'integer*4');
temp = fread(fid, j/8, 'real*8');
j = fread(fid, 1, 'integer*4');

%%save the matrix 
%%%%loop over the rows (temperature)
% filemark = 8*length(fr);
% for i = 1:m
%    fwrite(fid,filemark,'integer*4');
%    fwrite(fid,ks(i,:)','real*8');
%    fwrite(fid,filemark,'integer*4');
% end

kx = zeros(m,n);

for i = 1 : m
  j = fread(fid, 1, 'integer*4');
  kx(i,:) = fread(fid, [1,n], 'real*8');
  j = fread(fid, 1, 'integer*4');
end

fclose(fid);

freq = fstart : df : fend;
plot(freq,kx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc   = 0;
fr    = freq;
[m,n] = size(kx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('writing out the file .....')
fid = fopen(fname2, 'w', 'ieee-be');

fstart = fr(1);
fend   = fr(length(fr));
npts   = length(fr);
df     = fr(2)-fr(1);

%%header = f1,f2,df
filemark=8+8+8;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[fstart,fend,df],'real*8');
fwrite(fid,filemark,'integer*4');

%%header = CKD vers,m,n
filemark=4+4+4+4;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[loc,CKD,m,n],'integer*4');
fwrite(fid,filemark,'integer*4');

%save the temps
temp=100:10:400;
filemark = 8*length(temp);
fwrite(fid,filemark,'integer*4');
fwrite(fid,temp','real*8');
fwrite(fid,filemark,'integer*4');

%%save the matrix 
%%%%loop over the rows (temperature)
filemark = 8*length(fr);
for i = 1:m
   fwrite(fid,filemark,'integer*4');
   fwrite(fid,kx(i,:)','real*8');
   fwrite(fid,filemark,'integer*4');
   end
fclose(fid);

