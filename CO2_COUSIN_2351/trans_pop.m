%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trans_pop.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(pick_no);

disp('computing transition amplitudes and populations')
        a1=-0.2199475485E+01;   b1=0.9675055715E+00;
        c1=-0.8082711378E-03;   d1=0.2803987451E-05;
Q_296=a1+b1*296+c1*296^2+d1*296^3;
Q_t=a1+b1*temperature(i)+c1*temperature(i)^2+d1*temperature(i)^3;

trans_scale1=8*pi^3.*freqr.*(1-exp(-btz.*freqr...
	/296))/3/6.626176e-34/2.99792458e10.*9;
trans_scale2=(2.*jr+1).*exp(-btz.*elowerr/296).*1e-36/Q_296;
populationr=trans_scale1.*trans_scale2*1e-7;
transition_prob=strenr./populationr;  % at 296K
trans_amplr=sqrt(transition_prob);
trans_scale1_t=8*pi^3.*freqr.*(1-exp(-btz.*freqr...
	/temperature(i)))/3/6.626176e-34/2.99792458e10.*9;
trans_scale2_t=(2.*jr+1).*exp(-btz.*elowerr/temperature(i)).*1e-36/Q_t;
population_tr=trans_scale1_t.*trans_scale2_t*1e-7;

trans_scale1=8*pi^3.*freqp.*(1-exp(-btz.*freqp...
	/296))/3/6.626176e-34/2.99792458e10.*9;
trans_scale2=(2.*jp+1).*exp(-btz.*elowerp/296).*1e-36/Q_296;
populationp=trans_scale1.*trans_scale2*1e-7;
transition_prob=strenp./populationp;
trans_amplp=sqrt(transition_prob);
trans_scale1_t=8*pi^3.*freqp.*(1-exp(-btz.*freqp...
	/temperature(i)))/3/6.626176e-34/2.99792458e10.*9;
trans_scale2_t=(2.*jp+1).*exp(-btz.*elowerp/temperature(i)).*1e-36/Q_t;
population_tp=trans_scale1_t.*trans_scale2_t*1e-7;

clear a1 b1 c1 d1 Q_296 trans_scale1 trans_scale2 transition_prob
clear trans_scale1_t trans_scale2_t Q_t

end