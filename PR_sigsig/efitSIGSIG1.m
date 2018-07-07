function y=efit(bds,params)

jq      = params(:,1);
elowerq = params(:,2);

% This function will fit the lower state energy levels of the even 
% j states given in the HITRAN tape to the equation:
%
%	Erot= B*J(J+1) - D*[J(J+1)]^2 + Evib
%
% where Erot is the rotational energy and Evib is the vibrational
% energy.  Evib is given by:
%
%	Evib = v1(dv1) +v2(dv2) +v3(dv3)
%
% where vi is the vibration mode and dvi is the change in vi.
% 
% The program operates as follows:
% (1) Parameters are passed to efit.m
%	(a) evib is the starting value of Evib in the above equation
% 	(a) bs is the starting value of B in the above equation
%	(b) ds is the starting value of D in the above equation
%	(c) start2 is a matrix containing:
%		1) the even J values from the hittomat file and 
%		2) the corresponding lower state energies also
%		   obtained from the hittomat file.
%	(d) evib is the vibrational energy change for the band.
% (2) The difference between the energy calculated by the above
%     equation and the given energies is calculated.
% (3) This difference is minimized via a least squares fit by finding
%     the optimum values for B and D.

global quiet

evib = bds(1);   	
bs   = bds(2);	
ds   = bds(3);	
erot = (bs*jq.*(jq+1) -ds*(jq.*(jq+1)).^2);
y    = erot+evib-elowerq;

ind = find(jq<70);
y   = y(ind);

%%%if ((bs < 0) | (ds < 0))
if (bs < 0)
  y = 10000*ones(size(y));
  end

pflag = quiet;
%%% pflag = 1;
if pflag > 0
  bds
  clf
  subplot(211); plot(jq(ind),elowerq(ind)-evib,'*',jq(ind),erot(ind),'-')
	        ylabel('Energy');title('Fitting for Energies');
  subplot(212); plot(jq(ind),y)
	        xlabel('f');ylabel('Fitting Diff');
  pause(0.1)
  end


