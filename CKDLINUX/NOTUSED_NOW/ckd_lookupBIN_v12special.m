%script file to produce BINARY fortran look up table at 1.0 cm-1

%run this from the main SPECTRA directory
%this is to make HYBRID ckd2.4,0.0 as : call it CKD12 as it is a hybrid
%     data is best fit from 600-1575, 1625-3000 using CKD foreign 2.4
%     data is best fit from       1590-1610     using CKD foreign 0.0

%so from  500-3000, use self CKD2.4
%   from  500-1575, use for  CKD2.4
%   from 1625-3000, use for  CKD2.4
%   from 1575-1625, use a blend of foreign( CKD0 + CKD24)

path(path,'/home/sergio/SPECTRA');
 
%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%% note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do 
%%% last freq point, so it will stop at 3000 cm-1
loc=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
CKD=24;
[fr,ks]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,CKD,...
                           1,0,1,loc,'CKDLINUX/watercont');
[m,n]=size(ks);

CKD=12;     %sergio
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf12.bin';
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
CKD = 00;
[fr,kf00]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,...
			     CKD,0,1,1,loc,'CKDLINUX/watercont');
CKD = 24;
[fr,kf24]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,...
			     CKD,0,1,1,loc,'CKDLINUX/watercont');
c24 = ones(size(fr));

jjlt= find((fr >= 1575) & (fr < 1600));
slope = (0-1)/(1600-1575);
c24(jjlt) = slope*(fr(jjlt)-1575) + 1.0;

jjgt= find((fr >= 1600) & (fr < 1625));
slope = (1-0)/(1625-1600);
c24(jjgt) = slope*(fr(jjgt)-1600) + 0.0;

c00 = 1.0-c24;
plot(fr,c24)

%there are 31 temperatures in watercont file
c24 = ones(31,1)*c24;
c00 = ones(31,1)*c00;

kf = c00.*kf00 + c24.*kf24;

[m,n]=size(kf);

CKD=12;     %sergio
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor12.bin';
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

