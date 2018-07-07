%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% efit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=efit(bds,jr,elowerr)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor;
global trans_ampl population population_t t_rawdata;
global voi_back 
global K_back strenrt w_forr w_selfr

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

	evib=bds(1);   	
	bs=bds(2);	
	ds=bds(3);	
	erot=(bs*jr.*(jr+1) -ds*(jr.*(jr+1)).^2);
	y=erot+evib-elowerr;


pflag=0;
if pflag==1
	clg;hold off;
	subplot(211);plot(jr,elowerr-evib,'*',jr,erot,':')
	ylabel('Energy');title('Fitting for Energies');
	subplot(212);plot(jr,y)
	xlabel('f');ylabel('Diff');pause
end














