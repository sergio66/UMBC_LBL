%%script file to produce BINARY fortran look up table at 1.0 cm-1

%run this from the main SPECTRA directory
%this is to make HYBRID ckd2.4,Machado-Strow-Tobin
%     call it CKD 55 as it is a hybrid .... uses latest estimates from
%                    /WATER/CONTINUUM/MATFILES/MAY9_BESTFIT
%     data is best fit from 600-1300, 1800-3000 using CKD foreign 2.4
%     data is best fit from       1300-1800     using RAL data

path(path,'/home/sergio/SPECTRA');

disp('**********************************************************************');
disp('        doing hybrid of ckd 51, 55 : see airs_continuum_info.txt  ');
disp('**********************************************************************');

multipliers = load('airs_continuum_info.txt');

fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf51.bin';
[ks51, freq, temp] = contread(fname);
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor51.bin';
[kf51, freq, temp] = contread(fname);

fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf55.bin';
[ks55, freq, temp] = contread(fname);
fname = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor55.bin';
[kf55, freq, temp] = contread(fname);

airsfreq = multipliers(:,2);
whichone = multipliers(:,3);
factor   = multipliers(:,4);
mult51 = find(whichone == 51);
mult55 = find(whichone == 55);

for ii = 1 : length(freq)
  d51 = abs(airsfreq(mult51) - freq(ii)); 
  pp51 = find(d51 == min(d51)); d51 = min(d51); 

  d55 = abs(airsfreq(mult55) - freq(ii));  
  pp55 = find(d55 == min(d55)); d55 = min(d55); 

  if ((freq(ii) <= min(airsfreq)) | (freq(ii) >= max(airsfreq)))
    index51(ii) = 1;
    index55(ii) = 0;
    mult(ii) = 1.0;
  elseif (d51 < d55)
    index51(ii) = 1;
    index55(ii) = 0;
    mult(ii) = factor(mult51(pp51(1)));
  else
    index51(ii) = 0;
    index55(ii) = 1;
    mult(ii) = factor(mult55(pp55(1)));
    end
  end

for ii = 1 : 31
  ks(ii,:) = ks51(ii,:).*index51.*mult + ks55(ii,:).*index55.*mult;
  kf(ii,:) = kf51(ii,:).*index51.*mult + kf55(ii,:).*index55.*mult;
  end

fr = freq;
semilogy(freq,ks(20,:),freq,ks51(20,:),freq,ks55(20,:))
semilogy(freq,kf(20,:),freq,kf51(20,:),freq,kf55(20,:))

CKD=60;
%%% note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do  
%%% last freq point, so it will stop at 3000 cm-1 
loc=1; 

%%%%%%%%%%%%%%save the data tessa
%do self

[m,n]=size(ks);
semilogy(fr,ks); pause(0.1)

extension = [num2str(CKD) '.bin'];
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' extension];
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

fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf'  extension];
copier = ['!cp ' fname '  '  fname2];
eval([copier]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do foreign
[m,n]=size(kf);
semilogy(fr,kf); pause(0.1)

fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor'  extension];
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

fname2 = ['/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor'  extension];
copier = ['!cp ' fname '  '  fname2];
eval([copier]);

plot(fr,ks./ks51)