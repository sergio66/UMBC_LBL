addpath /home/sergio/SPECTRA/CKDLINUX
addpath /home/sergio/SPECTRA

%script file to produce BINARY fortran look up table at 1.0 cm-1
%run this from the main SPECTRA directory

path(path,'/home/sergio/SPECTRA/');

%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%%note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do last freq point, so it will stop
%%%at 3000 cm-1
CKD    = 32;
loc    = 0;
divide = 1;

clear topts;

topts.local  = loc;
topts.CKD    = CKD;
topts.ffin   = ffin;
topts.nbox   = nbox;
topts.divide = divide;

% [outwave,out_array]=run6watercontinuum(gasID,fmin,fmax,... 
%       ffin,fmed,fcor,fstep,... 
%       xnear,xmed,xfar,nbox,strength_far,strength_near,divide,CKD,... 
%       selfmult,formult,usetoth,local,profname); 
%%%% has been replaced by 
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{ffin,nbox,... 
%                divide,CKD,selfmult,formult,local}); 
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{topts}) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
%[fr,ks] = ...
%  run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,divide,CKD,...
%  1,0,1,loc,'../CKDLINUX/watercont');
topts.selfmult = 1.0;
topts.formult  = 0.0;
rmer = ['!/bin/rm CNTNM.OPTDPT WATER.COEF']; eval(rmer)
[fr,ks] = run8watercontinuum(1,100,3001,'../CKDLINUX/watercont',topts);
figure(1); semilogy(fr,ks); title('Self 32'); disp('ret'); pause

[m,n]=size(ks);

extension = [num2str(CKD) '.bin'];
fname  = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
fname  = ['/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf' extension];
fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf' extension];

fprintf(1,'opening %s \n %s \n',fname,fname2);
[exist(fname) exist(fname2)]
if exist(fname) ~= 0
  fprintf(1,'fname = %s for SELF already exists \n',fname)
  error('SELF outfile already exists')
end

%error('self32 WOW')
fid = fopen(fname,'w','ieee-le');

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
%[fr,kf] = ...
%   run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,divide,CKD,...
%   0,1,1,loc,'../CKDLINUX/watercont');
topts.selfmult = 0.0;
topts.formult  = 1.0;
rmer = ['!/bin/rm CNTNM.OPTDPT WATER.COEF']; eval(rmer)
[fr,kf] = run8watercontinuum(1,100,3001,'../CKDLINUX/watercont',topts);
figure(1); semilogy(fr,ks); title('Self 32'); 
figure(2); semilogy(fr,kf); title('Forn 32'); disp('ret'); pause

[m,n]=size(kf);

extension = [num2str(CKD) '.bin'];
fname  = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor' extension];
fname  = ['/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor' extension];
fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor' extension];
fprintf(1,'opening %s \n %s \n',fname,fname2);
[exist(fname) exist(fname2)]
if exist(fname) ~= 0
  fprintf(1,'fname = %s for FORN already exists \n',fname)
  error('FORN outfile already exists')
end  
%error('forn25 WOW')
fid = fopen(fname,'w','ieee-le');

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
