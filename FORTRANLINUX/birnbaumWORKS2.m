function [chi1]=x_birnbaum2(f,freq,w_tot,temperature,tau2)
%
%   x_birnbaum2(f,freq,w_tot,temperature,tau2)
%%
%%
%% this is same as BirnBaumWORKS, except that it abs abs(nu) to make things 
%% symmetric
%
%   Compute Birnbaum chi function
%

tau0=0.72/temperature;
dnu=f-freq;
zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2));
plot(zz);pause(0.1);
ex=exp(tau2*w_tot + tau0*abs(dnu));

%chi1=zz.*(-pi/2).*besselh_tobin(1,sqrt(-1)*zz).*ex;
chi1=zz.*(-pi/2).*real(besselh(1,1,sqrt(-1)*zz)).*ex;

%chi2=zz.*(-pi/2).*besselh_tobin(1,sqrt(-1)*zz);

%plot(f,chi1); 
%[freq w_tot tau2]
%pause(0.1);

