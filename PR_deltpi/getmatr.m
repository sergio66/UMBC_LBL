function [matr]=getmatr(band,prb,freqq)
%this function get the appropriate matrix for 
%us to use kfull/klor

matr=ones(4,8);
mn=floor(min(freqq)); mx=ceil(max(freqq));
%%%doVmixSimple uses abs(f-15) to see if it wants to do k=klor*0.5
%%%so i've used abs-16 to be safe
matr(:,1)=[mn-20.0 mn-16.0 mx+16.0 mx+20.0]';


