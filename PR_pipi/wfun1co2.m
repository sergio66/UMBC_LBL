function dif = wfun1co2(x_co2air,elowerr,wr_tape,jr,temperature,stuff)

btz=stuff.btz;
B0=stuff.B0;

band=stuff.band;

% This function works in conjunction with wfunco2.m to comput W
%
% Comments and details on computing W are found in wfunco2.m
%
% Variables are:
%  x_co2air:		a1,a2,a3
%
%    elowerr:		HITRAN energies
%
%  wr_tape:		wr are the temperature adjusted air broadened widths
%			for the R-branch lines. (Only even J's) [length(wr)]=25
%
% temperature:		temperature

a1=x_co2air(1); a2=x_co2air(2); a3=x_co2air(3);

no_lines=length(elowerr); 

energy_diff=ones(no_lines,1)*[elowerr]'-[elowerr]*ones(1,no_lines);
energy_diff=abs(energy_diff);

K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);

J=ones(no_lines,1)*(2*jr+1)'; 

K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

%now put zeros
if (band==2322)
  aa=2:2:length(jr); bb=1:2:length(jr);
  K(aa,bb)=0.0; K(bb,aa)=0.0;
  end

dif=K+diag(wr_tape);

