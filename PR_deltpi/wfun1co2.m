function dif=wfun1co2(x_co2air,elowerq,wq_tape,jq,temperature,stuff)

a1=x_co2air(1); a2=x_co2air(2); a3=x_co2air(3);

B0=stuff.B0;
btz=stuff.btz;

no_lines=length(elowerq);	

energy_diff=ones(no_lines,1)*[elowerq]'-[elowerq]*ones(1,no_lines);
energy_diff=abs(energy_diff);
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
J=ones(no_lines,1)*(2.*jq+1)';  % This is NOT J, it is (2J+1)
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

dif=K+diag(wq_tape);











