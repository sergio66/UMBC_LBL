%%%%%%%%%%%%%%%%%%%%%%%%%%  which_fit.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file detemines which parameters: beta, band strength 
%	multiplier, and constant frequency shift, to fit for
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
flag4(1)=input('bsm     :  ');

bstart(1)=input('beta    :  ');
bstart(2)=input('tau2    :  ');
bstart(3)=input('f_shift :  ');
bstart(4)=input('c0 (gauss)    :  ');
bstart(5)=input('cspan (gauss) :  ');
