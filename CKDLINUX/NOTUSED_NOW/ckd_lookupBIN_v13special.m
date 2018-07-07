%script file to produce BINARY fortran look up table at 1.0 cm-1

%run this from the main SPECTRA directory
%this is to make HYBRID ckd2.3TEMP : call it CKD13 as it is a hybrid
%     use self continuum CKD2.3 for entire database
%     use dave tobin's CKD2.3 foreign continuum ... obtained from the following
% -rw-r--r--    1 sergio   strow        1936 Mar 29 10:23 ckd2p3for.mat
% -rw-r--r--    1 sergio   strow       37237 Mar 29 10:20 frgn_2.4_0_2250.ps
% -rw-r--r--    1 sergio   strow      198912 Mar 29 10:19 cf0.fig

path(path,'/home/sergio/SPECTRA');
 
%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%% note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do 
%%% last freq point, so it will stop at 3000 cm-1
loc=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
CKD=23;
[fr,ks]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,CKD,...
                           1,0,1,loc,'CKDLINUX/watercont');
[m,n]=size(ks);

CKD=13;     %sergio ... this is special version
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf13.bin';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do foreign
CKD = 23;
[fr,kf23]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,...
                             CKD,0,1,1,loc,'CKDLINUX/watercont');
%there is no temperature dependance, so do
kf23 = kf23(1,:);

load ckd2p3for.mat
%bdatax starts at 830 cm-1, so augment the data
ii1 = find(fr < 830);
ii2 = find(fr > 2500);
ii = [ii1 ii2];
bdatax1 = [fr(ii)   bdatax];
bdatay1 = [kf23(ii) bdatay];

kftemp = interp1(bdatax1,bdatay1,fr);
kf = ones(31,1)*kftemp;

semilogy(bdatax1,bdatay1,fr,kf23,bdatax,bdatay,fr,kftemp)

[m,n]=size(kf);

CKD=13;     %sergio ... this is special version
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor13.bin';
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

