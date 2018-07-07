function [elower,jall] = efitter(jr,band,elowerr,elower,prb);  

% Creating the initial parameters:
%   jr are the hittomat j levels.
%   elowerr are the corresponding lower state energy levels.  Note that the 
%      lower state energies given on the tape are TOTAL energies, not just
%      rotational energies.
%   bds are the starting values for the fit.

global v12p1

%disp('fitting for missing odd energy levels')
bds = [630; 	% Evib
     0.39;      % B
     1.5e-7];	% D

bds = [0.02970517; 	% Evib
       0.35939840;      % B
       1.098497009];	% D

if ((prb == 'R')|(prb == 'r'))
  bds(1) = 1000;
end

if v12p1 > 0
  clear options; 
  options(1) = 0;
  ebd = leastsq('efitSIGSIG',bds,options,[],jr,elowerr);
else
  clear options; 
  options = optimset('MaxIter',1000);
  %ebd = lsqnonlin(@efitSIGSIG,bds,[],[],options,jr,elowerr)
  ebd = bds;
  params = [jr elowerr];
  ebd = fsolve(@(ebd) efitSIGSIG1(ebd,params),bds,options);
end

format long e

% Coefficients which minimize the difference between the functional form of the
% energy and the given values are returned in the variable ebd.
% Details of this step are shown in efit.m

% The energy levels for ALL J are calculated:

jall = [0:max(jr)+1]'; 

evib = ebd(1);	

%%% Now I will calculate the missing odd J levels and then stick in the 
%%% even tape values.  Again, notice that I will use the energies from the
%%% tape when possible.
elower =  (ebd(2)*jall.*(jall+1) - ebd(3)*(jall.*(jall+1)).^2);

%addpath /home/sergio/MATLABCODE
%keyboard_nowindow

elower = elowerr-evib;

%  elower will be used as input for calculating W
%  NoTe that elower is only the ROTATIONAL energy.
% elowerrf =  (ebd(2)*jr.*(jr+1) - ebd(3)*(jr.*(jr+1)).^2);

clear start2 bds
%save efittervar.mat elowerr elower jr
