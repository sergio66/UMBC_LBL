%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wfun1co2.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dif=wfun1co2(x_co2air,elowerq,wq_tape,jq,temperature,stuff)

B0=stuff.B0; 
btz=stuff.btz; 
beta=stuff.beta;

% This function works in conjunction with wfunco2.m to comput W
%
% Comments and details on computing W are found in wfunco2.m
%
% Variables are:
%  x_co2air:		a1,a2,a3
%
%    elowerq:		HITRAN energies
%
%     wq_tape:		wq are the temperature adjusted air broadened widths
%			for the Q-branch lines. (Only even J's) [length(wq)]=25
%
% temperature:		temperature

a1=x_co2air(1); a2=x_co2air(2); a3=x_co2air(3);

no_lines=length(elowerq); 

energy_diff=ones(no_lines,1)*[elowerq]'-[elowerq]*ones(1,no_lines);
energy_diff=abs(energy_diff);
  %
  % Energy gap scaling law 
  %
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
  %
  % Detailed balance
  %

J=ones(no_lines,1)*(2*jq+1)';    % This is NOT J, it is (2J+1)

K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

% The off diagonal elements of W are calculated by:
%
%  	Wj'j = beta*Kj'j 	 where Kj'j is for the pie state


%dif=beta*K+diag(wq_tape);      %%%%%%%%% originally in the Q-sig pi code
dif=K+diag(wq_tape);           %%%%%%%%% all the other codes have this
                               %%%%%%%%% eg Q-deltpi, PR-sigsig etc




