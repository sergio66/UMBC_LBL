%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wfun1co2.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dif=wfun1co2(x_co2air,elowerr,wq_tape,jr,num)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor;
global trans_ampl population population_t t_rawdata;
global voi_back
global K_back strenrt w_forr w_selfr

% This function works in conjunction with wfunco2.m to comput W
%
% Comments and details on computing W are found in wfunco2.m
%
% Variables are:
%  x_co2air:		a1,a2,a3
%
%    elowerr:		HITRAN energies
%
%  wq_tape:		wq are the temperature adjusted air broadened widths
%			for the Q-branch lines. (Only even J's) [length(wq)]=25
%
% temperature:		temperature

a1=x_co2air(1); a2=x_co2air(2); a3=x_co2air(3);

no_lines=length(elowerr); 

energy_diff=ones(no_lines,1)*[elowerr]'-[elowerr]*ones(1,no_lines);
energy_diff=abs(energy_diff);

K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature(num));

J=ones(no_lines,1)*(2*jr+1)'; 

K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature(num)).*J'./J;

dif=K+diag(wq_tape);

