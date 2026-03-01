%script file to produce BINARY fortran look up table at 1.0 cm-1
% copies the CKD1 data and applies Scott's multiplier

%run this from the main SPECTRA directory

path(path,'/home/sergio/SPECTRA/');

%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%%note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do last freq point, so it will stop
%%%at 3000 cm-1
CKD    = 4;
loc    = 0;
divide = 1;

clear topts;

FNAME = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.dat';
mult = load(FNAME);
fprintf(1,'mutiplier for CKD 4 = CKD 1 * mult is at \n %s \n',FNAME);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
ckd1 = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf1.bin';
[kx, fr, temp] = contread(ckd1);

multiplier = interp1(mult(:,2),mult(:,5),fr);
oo = find(fr <= mult(1,2) | fr >= mult(2378,2));
multiplier(oo) = 1.0;

multiplier = ones(31,1)*multiplier;

ks = kx.*multiplier;
[m,n]=size(ks);

extension = [num2str(CKD) '.bin'];
fname  = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
fid=fopen(fname,'w','ieee-le');

fstart = fr(1);
fend   = fr(length(fr));
npts   = length(fr);
df     = fr(2)-fr(1);

%%header = f1,f2,df
filemark=8+8+8;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[fstart,fend,df],'real*8');
fwrite(fid,filemark,'integer*4');

%%header = loc,CKD vers,m,n
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
   fwrite(fid,ks(i,:)','real*8');
   fwrite(fid,filemark,'integer*4');
   end
fclose(fid);

copier = ['!/bin/cp ' fname ' ' fname2];
fprintf(1,' the copy string is %s \n',copier);
eval(copier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do for
ckd1 = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor1.bin';
[kx, fr, temp] = contread(ckd1);

multiplier = interp1(mult(:,2),mult(:,5),fr);
oo = find(fr <= mult(1,2) | fr >= mult(2378,2));
multiplier(oo) = 1.0;

multiplier = ones(31,1)*multiplier;

kf = kx.*multiplier;
[m,n]=size(kf);

extension = [num2str(CKD) '.bin'];
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor' extension];
fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor' extension];
fid=fopen(fname,'w','ieee-le');

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
   fwrite(fid,kf(i,:)','real*8');
   fwrite(fid,filemark,'integer*4');
   end
fclose(fid);

copier = ['!/bin/cp ' fname ' ' fname2];
fprintf(1,' the copy string is %s \n',copier);
eval(copier);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor1.bin';
[kf1, fr, temp] = contread(ckd);
ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor2.bin';
[kf2, fr, temp] = contread(ckd);
ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor4.bin';
[kf4, fr, temp] = contread(ckd);

ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf1.bin';
[ks1, fr, temp] = contread(ckd);
ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf2.bin';
[ks2, fr, temp] = contread(ckd);
ckd = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf4.bin';
[ks4, fr, temp] = contread(ckd);

semilogy(fr,kf1); pause
plot(fr,kf4./kf1); pause
semilogy(fr,ks1); pause
plot(fr,ks4./ks1); pause

plot(fr,kf4./kf1,fr,kf4./kf2)

plot(fr,ks4./ks1,fr,ks4./ks2)
