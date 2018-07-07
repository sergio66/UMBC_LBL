%script file to produce BINARY fortran look up table at 1.0 cm-1

%run this from the main SPECTRA directory
%this is to make HYBRID ckd2.4,Machado-Strow-Tobin
%     call it CKD 52 as it is a hybrid .... uses latest estimates from
%                    /WATER/CONTINUUM/MATFILES/MAY13
%     data is best fit from 600-1300, 1800-3000 using CKD foreign 2.4
%     data is best fit from       1300-1800     using RAL data

path(path,'/home/sergio/SPECTRA');

%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%% note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do 
%%% last freq point, so it will stop at 3000 cm-1
loc=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
CKD=52;
quick.CKD  = 52;
quick.local= 0;
quick.ffin = ffin;
quick.nbox = nbox;
quick.divide   = 1; 
quick.selfmult = 1; 
quick.formult  = 0; 
[fr,ks]=run7watercontinuum(1,100,3001,'CKDLINUX/watercont',quick);
[m,n]=size(ks);
semilogy(fr,ks); pause

CKD=52;     %sergio
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf52.bin';
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

fname2 = '/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf52.bin';
copier = ['!cp ' fname '  '  fname2];
eval([copier]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do foreign
quick.CKD  = 52;
quick.local= 0;
quick.ffin = ffin;
quick.nbox = nbox;
quick.divide   = 1; 
quick.selfmult = 0; 
quick.formult  = 1; 
[fr,kf]=run7watercontinuum(1,100,3001,'CKDLINUX/watercont',quick);
[m,n]=size(kf);
semilogy(fr,kf); pause

CKD=52;     %sergio
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor52.bin';
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

fname2 = '/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor52.bin';
copier = ['!cp ' fname '  '  fname2];
eval([copier]);
