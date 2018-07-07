function [W_co2for,W_co2self] = ...
    wfunco2er(jr,elower,elowerr,w_selfr,w_forr,band,jall,temperature,stuff); 

%  This file does two things:
%	(1) Finds the coefficients a1, a2, & a3 of the energy gap
%	    scaling law by method of least-squares.
%	(2) Computes the K matrix and the W matrix.

global v12p1

% (1)
% Utilize "wfunco2.m" to find  a1,a2,a3:
% xx contains the starting values of the coefficients
xx = [ -6.202132765289518e-02   3.455577145107923e-01   1.009589099173803e+00];

%disp('fit for foriegn coefficients of energy-gap scaling law')
if v12p1 > 0
  clear options; 
  options(1) = 0; 
  options(5) = 0; options(7) = 1; options(14) = 60;
  a1a2a3_for  =  ...
    leastsq('wfunco2',xx,options,'wgradco2',elower,w_forr,jall,...
            temperature,stuff);
else
  clear options; 
  options = optimset('MaxIter',1000);
  a1a2a3_for = lsqnonlin(@wfunco2,xx,[],[],options,elower,w_forr,jall,...
                       temperature,stuff);
  end
W_co2for = wfun1co2(a1a2a3_for,elowerr,w_forr,jr,temperature,stuff);

% Now do the same for self-broadened widths.
%disp('fit for self coefficients of energy-gap scaling law')
if v12p1 > 0
  a1a2a3_self = leastsq('wfunco2',xx,options,'wgradco2',...
              elower,w_selfr,jall,temperature,stuff);
else
  a1a2a3_self = lsqnonlin(@wfunco2,xx,[],[],options,...
                        elower,w_selfr,jall,temperature,stuff);
  end
W_co2self = wfun1co2(a1a2a3_self,elowerr,w_selfr,jr,temperature,stuff);

clear options














