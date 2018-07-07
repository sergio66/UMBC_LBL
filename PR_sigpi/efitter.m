function [elower,jall] = efitter(jq,band,elowerq,elower,prb); 

%%%%%%%%%%%%%%%%%%%%%%%%%  efitter.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program works in conjunction with loader.m and efit.m as follows:
% (1) Hittomat data is loaded and the pertinent Q-branch lines are selected
%     out by "loader.m"
% (2)"efitter.m" is called to calculate the lower state rotational energies
%	This is necessary because these values are needed to calculate the 
%	state-to-state collision rates.  The odd j energy levels are not
%	in several Q-branches because of symmetry requirements.  (These 
%	levels are, however, present in the P and R branches of the same
%	vibratonal transitions for some of these Q-branches.)  To standardize
%	the solution of finding these odd energy levels, a least squares
%	fit will be used to find a functional form for the given even energy
%	levels and then the missing odd levels can be calculated.  Although
%	a more precise way to get these levels would be to use the odd level
%	values given in the P & R branch bands, the least-squares fit is more
%	standard because the odd levels are not present in all cases.  In
%	addition, the precision required for the values of these energy levels
%	is well met by the least-squares method.
%
% Creating the initial parameters:
%   jq are the hittomat even j levels. (jq =  2,4,6,8,...50)
%   elowerq are the corresponding lower state energy levels.  Note that the 
%      lower state energies given on the tape are TOTAL energies, not just
%      rotational energies.
%   bds are the starting values for the fit.

global quiet v12p1

if quiet > 0
  fprintf(1,'fitting for missing odd energy levels \n')
  end

bds = [670; 	% Evib
     0.39;		% B
     1.5e-7];	% D

if v12p1 > 0
  clear options; 
  options(1) = 0;
  ebd = leastsq('efit',bds,options,[],jq,elowerq);
else
  clear options; 
  options = optimset('MaxIter',1000);
  %ebd = lsqnonlin('efit',bds,[],[],options,jq,elowerq);
  ebd = bds; 
  params = [jq elowerq]; 
  ebd = fsolve(@(ebd) efit1(ebd,params),bds,options); 
  end

% Coefficients which minimize the difference between the functional form of the
% energy and the given values are returned in the variable ebd.
% Details of this step are shown in efit.m

% The energy levels for ALL J (except J=0) are calculated:
if ((band  ==  720) | (band  ==  791))
  if (prb == 'R')
    if (rem(max(jq),2)==0)
      jall = [0:max(jq)]'; 
    else 
      jall = [0:max(jq)+1]'; 
      end
  elseif  (prb == 'P')
    if (rem(max(jq),2)==0)
      jall = [0:max(jq)]'; 
    else 
      jall = [0:max(jq)+1]'; 
      end
  elseif (prb == 'Q')
  %this is directly from Q branch code; use also for P and R branches
    jall = [0:max(jq)]';  % 0 to 50
    end

elseif band  ==  667
  if (prb == 'R')
    if (rem(max(jq),2)==0)
      jall = [0:max(jq)+1]'; 
    else 
      jall = [0:max(jq)]'; 
      end
  elseif  (prb == 'P')
    if (rem(max(jq),2)==0)
      jall = [0:max(jq)]'; 
    else 
      jall = [0:max(jq)+1]'; 
      end
  elseif (prb == 'Q')
  %this is directly from Q branch code; use also for P and R branches
    jall = [0:max(jq)]';  % 0 to 50
    end

end

evib = ebd(1);	

%%% Now I will calculate the missing odd J levels and then stick in the 
%%% even tape values.  Again, notice that I will use the energies from the
%%% tape for the even levels.
elower =  (ebd(2)*jall.*(jall+1) - ebd(3)*(jall.*(jall+1)).^2);
index_e = jq+1;
elower(index_e) = elowerq-evib;

%  elower will be used as input for calculating W
%  *******NoTe that elower is only the ROTATIONAL energy.*******

eloweref =  (ebd(2)*jq.*(jq+1) - ebd(3)*(jq.*(jq+1)).^2);

clear start2 bds eloweref
