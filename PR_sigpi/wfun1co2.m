%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wfun1co2.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dif=wfun1co2(x_co2air,elowerq,wq_tape,jq,temperature,stuff)

B0=stuff.B0; 
btz=stuff.btz; 
beta=stuff.beta;

%%%% note that this is directly the code from the Q branch
%%%% this is possible because wfunco2.m does the fitting of a1 a2 a3 to
%%%% the linewidths, usimg jall ie using the coupling between all the
%%%% lower and upper lines j=0,1,2,3,4,5, ...
%%%% while this file uses jq=0,2,4,6, etc

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


%%%%%%dif=beta*K+diag(wq_tape);  %%%% not doing this
dif=K+diag(wq_tape);
%%%%%%%%% this is to test direct comparison to Lorentz 
%%%%%%%%  dif=diag(wq_tape);



