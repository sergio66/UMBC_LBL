%%%%%%%%%%%%%%%%%%%%%%%%% wfunco2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dif=wfun(xx,elower,wq_tape,jall,prb,num)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor;
global trans_ampl population population_t t_rawdata;
global voi_back
global K_back strenrt w_forr w_selfr

% Variables are:
%  	 xx:		starting values for a1,a2,a3
%
%    elower:		HITRAN energies
%
%     wq_tape:		wq are the temperature adjusted widths
%			for the R-branch lines. (Only even J's)
%
% temperature:		temperature

no_lines=length(elower);  

%%%%%%%%%%%%%%%%%%%%%% Energy difference calculation %%%%%%%%%%%%%%%%%%%%%%%
% Here, the energy difference between levels is calculated

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff);

%%%%%%%%%%%%%%%%%%%%%% Calculation of collision rates %%%%%%%%%%%%%%%%%%%%
a1=xx(1);a2=xx(2);a3=xx(3);
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature(num));

  % Detailed balance is achieved by using the equation:
  %   Kjj' *(2j'+1)*exp(-Ej'/kT) = Kj'j *(2j+1)*exp(-Ej/kT)
J=ones(no_lines,1)*(2*jall+1)';    % This is NOT J, it is 2J+1

K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature(num)).*J'./J;

%%%%%%%%% Calculating the widths: (i.e. the diagonal elements of W) %%%%%%%%
%
%          index    1   2   3   4
%   K=           ___________________
%            1   | k00 k10 k20 k30 | 0
%	     2	 | k01 k11 k21 k31 | 1
%	     3	 | k02 k12 k22 k32 | 2
%	     4	 | k03 k13 k23 k33 | 3
%		 ___________________
%                   0   1   2   3    J

i_even=find(rem(jall,2)==0);             % i=1,3,5...   J=0,2,4...
i_odd=find(rem(jall,2)~=0);              % i=2,4,6...   J=1,3,5...

	width_even_lower=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
	width_odd_upper=-0.5*sum(K(i_odd,i_odd));       % k11+k13+k15...
if prb=='P'
	width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
end

widthq=width_even_lower + width_odd_upper;

%%% Calculated the difference between fitted widths and data widths%%%%%%%%%%
dif=(widthq'-wq_tape)';

pflag=0;
if pflag==1
	clg;
	jqplot=(0:2:max(jall));
	if prb=='P'
		jqplot=(2:2:max(jall));
	end
	subplot(2,1,1);plot(jqplot,wq_tape,'*',jqplot,widthq')
	ylabel('Width');title('Fitting for a1,a2,a3 via widths');
	subplot(2,1,2);plot(jqplot,dif)
	xlabel('f');ylabel('Diff.');drawnow
end




