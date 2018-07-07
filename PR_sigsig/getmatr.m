function [matr]=getmatr(band,prb,freqq)
%this function get the appropriate matrix for us to use kfull/klor

%done from findratio.m

%this is the strongest band, so do everything CORRECT!!!!!!!!
matr=zeros(4,8);
matr(1,1)=min(freqq)-110;  
matr(2,1)=min(freqq)-120;
matr(3,1)=max(freqq)+110;  
matr(4,1)=max(freqq)+120;
