function [Y_1st]=y1s(jq,freqq,elowerq,strenq,W_matrix,trans_ampl,stuff,beta)

%   This program computes the first order mixing coefficients: Y_1st.
%   
%   They are calculated with the equation:
%
%	Yi= 2* SUM{(dk/di)*Wki/(vi-vk)} 
%		where the sum is over k not equal to i
%		with dk,i= dipole matrix elements for lines k,i
%		     Wki = W matrix element for relaxation from i to k
%		     vk,i= wavenumber of line k,i
%

%no_linesq=length(freqq);  
%Y_1stmix=zeros(1,no_linesq);
% for n=1:no_linesq
%  for m=1:no_linesq
%   if m~=n
%    Y_1stmix(n)=Y_1stmix(n)+trans_ampl(m)*(W_matrix(m,n))/(freqq(n)-freqq(m));
%   end
%  end
% end
% Y_1st=2*Y_1stmix./trans_ampl';

% Note that in this sum, trans_ampl(n) is not summed over so we can pull it
% out of the sum.

IncludeMaxer
[m,n]=size(W_matrix);
if (m ~= n) 
  error('need a square Wmatrix!!!'); 
  end 
W_matrix0=zeros(MaxPQR,MaxPQR); 
W_matrix0(1:m,1:m)=W_matrix; 
Y_1stmixN=doFindMix(trans_ampl,W_matrix0,freqq);
Y_1st=Y_1stmixN';

%%%%% wanna compare to ordering stuff
IncludeMaxer
[m,n]=size(W_matrix);
if (m ~= n) 
  error('need a square Wmatrix!!!'); 
  end 
freqq0        = freqq;
trans_ampl0   = trans_ampl;
[freqq,ii]=sort(freqq);
trans_ampl=trans_ampl(ii);
B=W_matrix;
for jj=1:length(ii)
  W_matrix(:,jj)=B(:,ii(jj));
  end
W_matrix0=zeros(MaxPQR,MaxPQR); 
W_matrix0(1:m,1:m)=W_matrix; 
Y_1stmixN=doFindMix(trans_ampl,W_matrix0,freqq);
Y_1st=Y_1stmixN';
freqq        = freqq0;
trans_ampl   = trans_ampl0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this part is slghtly different from y1s_WORKS_july12_2015.m as we need to re-generate first order linemix tables

iNLTE_Print = -1;
%% need to be correct here, consistent with topts.hitran <<<<<<<<<<<<<<<<<
%% need to be correct here, consistent with topts.hitran <<<<<<<<<<<<<<<<<
%% need to be correct here, consistent with topts.hitran <<<<<<<<<<<<<<<<<
HITRAN      = 2008;  co2ppm = 370;            
HITRAN      = 2000;  co2ppm = 370;
HITRAN      = 2000;  co2ppm = 385;
HITRAN      = 2012;  co2ppm = 385;

if iNLTE_Print > 0
  disp(' ')
  disp(' >>>> !!! >>>> save first order linemix coeffs : HOPE you set HITRAN and co2ppm in CO2_COMMON/y1s.m correctly!!!')
  disp(' >>>> !!! >>>> save first order linemix coeffs : HOPE you set HITRAN and co2ppm in CO2_COMMON/y1s.m correctly!!!')
  disp(' >>>> !!! >>>> save first order linemix coeffs : HOPE you set HITRAN and co2ppm in CO2_COMMON/y1s.m correctly!!!')
  disp(' ')
  
  %% see /home/sergio//KCARTA/SRC/NONLTE/M_Files_for_kcarta_NLTE_LBL_runs/USUALLAYERS/drivelineparameters.m
  %% new  post-July 2015
  %% see n_gas_wt_spectra.f, subr NLTEBandMapper
  %% strong bands
  stoptsRMV.iKCBand = [2310 2311     2320 2321 2322   2350 2351 2352 2355   2353 2253     2354 2254];
  stoptsRMV.iISO    = [   1    2       1    2    3      1    2    3    4      1   2        1    2  ];
  stoptsRMV.iL      = [   4    4       2    2    2      1    1    1    1      3   3        5    5  ];
  stoptsRMV.iU      = [  24   24      16   16   16      9    9    9    9      23  23       25  25  ];

  %other not so strong bands, but also in NLTE
  wtoptsRMV.iKCBand = [2110 2120 2130 2140 2150 2160 2170 2180];
  wtoptsRMV.iISO    = [  1   1    1    1    1    1    1    1  ];
  wtoptsRMV.iL      = [  2   3    4    5    3    6    7    8  ];
  wtoptsRMV.iU      = [ 15  25   22   23   22   36   37   38  ];

  toptsRMV.iKCBand = [stoptsRMV.iKCBand wtoptsRMV.iKCBand];
  toptsRMV.iISO    = [stoptsRMV.iISO    wtoptsRMV.iISO   ];
  toptsRMV.iL      = [stoptsRMV.iL      wtoptsRMV.iL     ];
  toptsRMV.iU      = [stoptsRMV.iU      wtoptsRMV.iU     ];

  outfile = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX_TEST/';
  if HITRAN ~= 2000
    hitx = num2str(HITRAN);
    hitx = hitx(3:4);
    hitstr = num2str(hitx);
  else
    hitstr = '2k';
  end
  outfile = ['/asl/data/kcarta_sergio/KCDATA/NLTE/AIRSCO2/LINEMIX_H' hitstr '_' num2str(co2ppm) '/PQR/'];
  
  outfile = [outfile num2str(round(stuff.temperature)) 'K_band' num2str(stuff.band)];
  outfile = [outfile  '_linemix' stuff.prb '.dat'];
  %outfile = [outfile  '_linemix.dat'];  

  fprintf(1,'band = %4i, T = %8.5f : saving first order limemix coeffs to  = %s \n',...
             stuff.band,round(stuff.temperature),outfile);
  if isfield(stuff,'beta')	     
    fprintf(1,'beta, doc = %8.6f %8.6f \n',stuff.beta,stuff.duration)
  else
    disp('no beta for this band !!! ???')
  end
  
  plot(freqq,Y_1st,'o-'); grid; xlabel('freqs'); ylabel('Y1 = first order linemix');
  title(['T = ' num2str(stuff.temperature) 'K, band = ' num2str(stuff.band)])
  pause(0.1)
  
  [haha,ii] = sort(freqq);
  fafa=freqq(ii)';
  yaya=Y_1st(ii);
  %popo = [1:length(freqq); fafa; Y_1st(ii)];
  popo = [jq'; fafa; Y_1st(ii)];  
  if (fafa(1) > fafa(length(fafa)))
    popo = fliplr(popo);
  end

  %fprintf(1,'%4i %8.6e   %8.6e \n',popo);

  gasID = 2;
  iNum  = length(freqq);
  iXYZ = find(toptsRMV.iKCBand == stuff.band);
  if length(iXYZ) == 1
    iISO  = toptsRMV.iISO(iXYZ);
    jLow  = toptsRMV.iL(iXYZ);
    jHigh = toptsRMV.iU(iXYZ);
  else
    iISO = -1;
    jLow = -1;
    jHigh = -1;
  end
  fid = fopen(outfile,'w');
  fprintf(fid,'%2i %3i %2i %2i %2i \n',gasID,iNum,iISO,jLow,jHigh);
  fprintf(fid,'%4i %8.6e   %8.6e \n',popo);
  fclose(fid);
  
end

%%%%also note that the code in pipi/WORKS, deltdelt/WORKS, for some strange 
%%%%reason, treated the P branch different :
%function Y_1st=y1s(jr,freqr,elowerr,strenr,W_matrix,trans_ampl,stuff,beta,prb)
%if prb=='P'
%  temp=zeros(length(Y_1st),1);
%  for i=1:length(Y_1st)
%    temp(i)=Y_1st(length(Y_1st)-i+1);
%    end
%  Y_1st=temp;
%  end
%this flipping caused trouble : kmix/klor >> 1 or <<<< 0 !!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






