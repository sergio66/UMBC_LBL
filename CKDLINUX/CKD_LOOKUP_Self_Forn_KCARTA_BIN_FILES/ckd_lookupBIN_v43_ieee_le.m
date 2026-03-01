addpath /home/sergio/SPECTRA
addpath /home/sergio/SPECTRA/CKDLINUX

%script file to produce BINARY fortran look up table at 1.0 cm-1
%run this from the main SPECTRA directory

path(path,'/home/sergio/SPECTRA/');
path(path,'/home/sergio/SPECTRA/CKDLINUX/');

%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin = 0.0025;
nbox = 1/ffin;

%%%note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do last freq point, so it will stop
%%%at 3000 cm-1
CKD    = 43;
loc    = 0;
divide = 1;

clear topts;

topts.local  = loc;
topts.CKD    = CKD;
topts.ffin   = ffin;
topts.nbox   = nbox;
topts.divide = divide;

extension = [num2str(CKD) '.bin'];
profile_file_for_run8 = '/home/sergio/SPECTRA/CKDLINUX/watercont';

% [outwave,out_array]=run6watercontinuum(gasID,fmin,fmax,... 
%       ffin,fmed,fcor,fstep,... 
%       xnear,xmed,xfar,nbox,strength_far,strength_near,divide,CKD,... 
%       selfmult,formult,usetoth,local,profname);
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [fr,kf] = ...
%   run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,divide,CKD,...
%   0,1,1,loc,profile_file_for_run8);
%
%[fr,ks] = ...
%  run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,divide,CKD,...
%  1,0,1,loc,profile_file_for_run8);
%
%%%% has been replaced by
%
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{ffin,nbox,... 
%                divide,CKD,selfmult,formult,local}); 
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{topts}) 
%
% [fr,k] = run8watercontinuum(gasID,fmin,fmax,profname,{ffin,nbox,... 
%                divide,CKD,selfmult,formult,local}); 
% [fr,k] = run8watercontinuum(gasID,fmin,fmax,profname,{topts}) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SELF 
topts.selfmult = 1.0;
topts.formult  = 0.0;
rmer = ['!/bin/rm CNTNM.OPTDPT WATER.COEF']; eval(rmer)
[fr,ks] = run8watercontinuum(1,100,3001,profile_file_for_run8,topts);
figure(1); semilogy(fr,ks); title('Self 43'); disp('ret'); pause

[m,n] = size(ks);

%% fname      = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
%  fname      = ['/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf' extension];
%  fname2copy = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf' extension];

fname      = ['/home/sergio/asl/rta/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
fname2copy = ['/home/sergio/asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf' extension];

fprintf(1,'opening %s \n %s \n',fname,fname2copy);
[exist(fname) exist(fname2copy)]
if exist(fname) ~= 0
  fprintf(1,'fname = %s for SELF already exists \n',fname)
  error('SELF outfile already exists')
end
if exist(fname2copy) ~= 0
  fprintf(1,'fname2copy = %s for SELF already exists \n',fname2copy)
  error('SELF outfile2copy already exists')
end

%error('self43 WOW')

kwriteout = ks;
fwrite_to_CKD_S_or_F_bin

copier = ['!/bin/cp ' fname ' ' fname2copy];
fprintf(1,' the copy string is %s \n',copier);
eval(copier);

lser = ['!ls -lt ' fname ' ' fname2copy];
eval(lser)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FOREIGN
topts.selfmult = 0.0;
topts.formult  = 1.0;
rmer = ['!/bin/rm CNTNM.OPTDPT WATER.COEF']; eval(rmer)
[fr,kf] = run8watercontinuum(1,100,3001,profile_file_for_run8,topts);
figure(1); semilogy(fr,ks); title('Self 43'); 
figure(2); semilogy(fr,kf); title('Forn 43'); disp('ret'); pause

[m,n] = size(kf);

%% fname  = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor' extension];
%  fname  = ['/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor' extension];
%  fname2copy = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor' extension];

fname      = ['/home/sergio/asl/rta/kcarta/KCARTADATA/General/CKDieee_le/CKDFor' extension];
fname2copy = ['/home/sergio/asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor' extension];

fprintf(1,'opening %s \n %s \n',fname,fname2copy);
[exist(fname) exist(fname2copy)]
if exist(fname) ~= 0
  fprintf(1,'fname = %s for FORN already exists \n',fname)
  error('FORN outfile already exists')
end  
if exist(fname2copy) ~= 0
  fprintf(1,'fname2copy = %s for FORN already exists \n',fname2copy)
  error('FORN outfile2copy already exists')
end

%error('forn43 WOW')

kwriteout = kf;
fwrite_to_CKD_S_or_F_bin

copier = ['!/bin/cp ' fname ' ' fname2copy];
fprintf(1,' the copy string is %s \n',copier);
eval(copier);

lser = ['!ls -lt ' fname ' ' fname2copy];
eval(lser)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
