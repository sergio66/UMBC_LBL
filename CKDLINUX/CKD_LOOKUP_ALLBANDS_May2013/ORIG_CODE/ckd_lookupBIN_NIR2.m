% script file to produce BINARY fortran look up table at 1.0 * FACTOR cm-1
% run this from the main SPECTRA directory

%% technically, CKD1 == CKD6 for 140-310 cm-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% aalways want 2901 points for kCARTA to read (gotta fix this bug)
nptsWANT = 2901;   %%% from fmin to fmax + df eg 100-3000 at 1cm-1
fcenter = 6000; fmin  = 3000; fmax  =  8800; kaTag = 'm';  
deltav = 2.0;
foutdir = '/asl/data/kcarta/H2004_otherbands.ieee-le/NIR3550_5550/etc.ieee-le/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% eqn is ffin = (fmax+ffin-fmin)/2901 ==> 2900 ffin = (fmax-fmin)
%%%                                     ==> ffin = (fmax-fmin)/2900
ffin = (fmax+deltav-fmin)/2901

nbox = 1;

%% note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do 
%% last freq point, so it will stop at 3000 cm-1

CKD = 1;
CKD = input('Enter CKD : ');
loc = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theprofname = 'IPFILES/watercont';

%do self
quick.CKD  = CKD; 
quick.local= loc; 
quick.ffin = ffin; 
quick.nbox = nbox; 
quick.divide   = 1;  
quick.selfmult = 1;  
quick.formult  = 0;  
[fr,ks]=run8watercontinuum(1,fmin,fmax+ffin,theprofname,quick); 

[m,n]=size(ks);

fout_s = [foutdir 'CKDSelf' num2str(CKD) '.bin'];
fid = fopen(fout_s,'w','ieee-le');

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

clear quick 
%do forn 
quick.CKD  = CKD; 
quick.local= loc; 
quick.ffin = ffin; 
quick.nbox = nbox; 
quick.divide   = 1;  
quick.selfmult = 0;  
quick.formult  = 1;  
[fr,kf]=run8watercontinuum(1,fmin,fmax+ffin,theprofname,quick); 

[m,n]=size(kf);

fout_s = [foutdir 'CKDFor' num2str(CKD) '.bin'];
fid = fopen(fout_s,'w','ieee-le');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; 

iLog = input('normal (-1) or semilogy (+1) plot? ');
if iLog < 0
  h11 = subplot(211); plot(fr,ks);  ylabel('ks');
  h12 = subplot(212); plot(fr,kf);  ylabel('kf');
else
  h11 = subplot(211); semilogy(fr,ks);  ylabel('ks');
  h12 = subplot(212); semilogy(fr,kf);  ylabel('kf');
  end
adjust21(h11,h12,'even'); xlabel('wavenumber cm-1')


fprintf(1,'fmin,fmax = %6i %6i \n',fmin,fmax);
[m,n1]=size(ks); fprintf(1,'have %5i x %5i for ks \n',m,n1);
[m,n2]=size(kf); fprintf(1,'have %5i x %5i for kf \n',m,n2);
if n1 ~= 2901 | n2 ~= 2901
  disp(' oh oh kcarta always needs 2901 points ...') 
  end