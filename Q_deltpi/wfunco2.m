function dif=wfun(xx,elower,wq_tape,jall,beta_pi,beta_delt,temperature,stuff)

global quiet

no_lines=length(elower);	

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff);

B0=stuff.B0;
btz=stuff.btz;

a1=xx(1);a2=xx(2);a3=xx(3);
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
J=ones(no_lines,1)*(2.*jall+1)';    % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

i_even=find(rem(jall,2)==0);   % i=2,4,6...   J=2,4,6...
i_odd=find(rem(jall,2)~=0);    % i=1,3,5...   J=1,3,5...
i_even1=i_even-1;              
i_odd2=i_even;
i_odd3=i_odd(2:length(i_odd)); % i=3,5,7...   J=3,5,7...

width_even_lower(i_even1)=-0.5*beta_delt*sum(K(i_even,i_even));
width_even_lower(i_odd2)=-0.5*(1-beta_delt)*sum(K(i_even,i_odd3));
width_odd_lower(i_even1)=-0.5*(1-beta_delt)*sum(K(i_odd3,i_even));
width_odd_lower(i_odd2)=-0.5*beta_delt*sum(K(i_odd3,i_odd3));

width_even_upper(i_even1)=-0.5*beta_pi*sum(K(i_even,i_even));
width_even_upper(i_odd2)=-0.5*(1-beta_pi)*sum(K(i_even,i_odd3));
width_odd_upper(i_even1)=-0.5*(1-beta_pi)*sum(K(i_odd,i_even));
width_odd_upper(i_odd2)=-0.5*beta_pi*sum(K(i_odd,i_odd3));

widthq=width_even_upper+width_even_lower+width_odd_upper+width_odd_lower;

widthq=widthq(1:length(wq_tape));
dif=(widthq'-wq_tape)';

pflag=quiet;
if (pflag > 0)
  lq=1:length(wq_tape);
  plot(lq,widthq,lq,wq_tape); pause(0.05);
  end

pflag=quiet;
if pflag > 0
	%plot the results
	clg;
	jqplot=(2:max(jall));
	subplot(211);plot(jqplot,wq_tape,'*',jqplot,widthq)
	ylabel('Width');title('Fitting for a1,a2,a3 via widths');
	subplot(212);plot(jqplot,-dif)
	xlabel('f');ylabel('Diff.')
end

