function [Y_1st]=y1s(jq,freqq,elowerq,strenq,W_matrix,trans_ampl,stuff,beta)

%%% this MEX-calling version is copied from /SPECTRA/CO2_COMMON/y1s.m
%%% as need a square matrix, of size MaxLEn,Maxlen for the F77 mex file to work

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

%%%%%wanna compare to ordering stuff .... for the fits, DO NOT DO THIS
iorder = -1;
if iorder > 0
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






