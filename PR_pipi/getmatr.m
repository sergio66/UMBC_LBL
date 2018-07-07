function [matr]=getmatr(band,prb,freqq)
%this function get the appropriate matrix for us to use kfull/klor

matr=ones(4,8)*1e-10;        %ratios are typically 3e-2
matr(1,1)=min(freqq)-60;  
matr(2,1)=min(freqq)-70;
matr(3,1)=max(freqq)+160;  
matr(4,1)=max(freqq)+170;

