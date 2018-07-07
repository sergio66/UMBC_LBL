function [trans_ampl,population_t]=... 
              trans_pop(temperature,freqr,jr,elowerr,strenr,stuff) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trans_pop.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This program computes the transition amplitudes and populations
%   This calculation comes fromthe AGLF paper.
%
%   The basic formula to know is:
%	    line-strength(T)=population(T)*transition_amplitude
%   Transition amplitude are independent of temperature and therefore
%   it doesn't matter what temperature you calculate them at.  By convention
%   they will be calculated here at T=296.  Note that the strengths given on
%   the tape are for 296 and therefore no strength adjustment is necessary
%   here.   When calculating 1st-order mixing, the quantity d(k)/d(i), where
%   the d's are transition amplitudes, is needed.  Also note that from the
%   equation in the AGLF paper you can get the relation between strengths at 
%   different temperatures.
%
%   The population(T), population_t, will however, be used in the calculation
%   of full-mixing and therefore must be calculated at the correct temperature.
%
%   NOTE that transition amplitudes are computed at T=296K !!!
%   These will be used for calculating Y's and full-mixing absorption coef.
%
%   population_t is the population at the desired temperature and
%   will be used in calculating full mixing.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%there are differences between the original code and what the HITRAN manual has
% -- orig code had a mystery factor of 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

   % These are the coefficients for calculationg the total internal
   % partition sums.  They are found in /parsnip/scratch/Genln2/
   % Genln2/qtips.f  They are for CO2 (iso=626) 70-400K range

isotope         = stuff.isotope;

% NOTE that the stimulated emiision terms ARE included.
% These are the partition function constants for CO2 obtained from qtips.pas

[a1,b1,c1,d1]   = qqttiippss(isotope);

B0              = stuff.B0;  
btz             = stuff.btz;  
temperature_ref = stuff.temperature_ref;

%%%%%%%% Calculate transition amplitudes at T=296.
%%%%%%%% Therefore use strengths at 296 ( from the tape.)

kappa = (8*pi^3)/(3*6.626176e-34*2.99792458e10);
T     = temperature;

%%%%%%%%%% there is a mystery factor of "9" that we are removing from code
%trans_scale1=8*pi^3.*freqr.*(1-exp(-btz.*freqr...
%	/296))/3/6.626176e-34/2.99792458e10.*9;
T               = 296;         
Q               = a1+b1*T+c1*T^2+d1*T^3;
trans_scale1    = kappa*freqr.*(1-exp(-btz.*freqr/T));
trans_scale2    = (2.*jr+1).*exp(-btz.*elowerr/T).*1e-36/Q;
population      = trans_scale1.*trans_scale2*1e-7;
transition_prob = strenr./population;
trans_ampl      = sqrt(transition_prob);

%%%%% Calculate populations at desired temperature

%%%%%%%%%% there is a mystery factor of "9" that we are removing from code
%trans_scale1_t=8*pi^3.*freqr.*(1-exp(-btz.*freqr...
%	/temperature))/3/6.626176e-34/2.99792458e10.*9;
T              = temperature; 
Q              = a1+b1*T+c1*T^2+d1*T^3;
trans_scale1_t = kappa*freqr.*(1-exp(-btz.*freqr/T));
trans_scale2_t = (2.*jr+1).*exp(-btz.*elowerr/T).*1e-36/Q;
population_t   = trans_scale1_t.*trans_scale2_t*1e-7;

%%%% NoTe that if these populations are used for full mixing they should
%%%% really be divided by freqr to remove the stimulated emission terms. 
clear a1 b1 c1 d1 trans_scale1 trans_scale2 transition_prob
clear trans_scale1_t trans_scale2_t Q
