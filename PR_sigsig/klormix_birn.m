function [lor,lormix,lor_birn,lormix_birn]=...
          klormix(f,freqr,W_plus,jr,temperature,beta,... 
      w_selfr,w_forr,trans_ampl,population_t,stuff,strenrt,ymix);  
 
frequency_shift=stuff.frequency_shift; 
density=stuff.density; 
pressure_self=stuff.pressure_self; 
pressure_ref=stuff.pressure_ref; 
temperature_ref=stuff.temperature_ref; 
path_length=stuff.path_length; 
bsm=stuff.bsm; 
pressure_for=stuff.pressure_for; 
duration=stuff.duration; 

 %%%%%%%%%%%%%%%%%%%%%%%%%%% klormix_birn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes:
%		(1) Lorentz with 1st order mixing absorption coefficient
%			with and without Birnbaum Chi factor
%		(2) Lorentz absorption coefficient
%			with and without Birnbaum Chi factor
%	It also checks SUM{Si*Yi} to see if detailed balance worked.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('computing absorption coeffs using first order theory.')
freqr_shift=freqr+frequency_shift/100;
K_scale_lor=density*pressure_self/pressure_ref*temperature_ref*path_length...
		/temperature/pi;

lor=zeros(length(f),1);		  % lorentz
lor_birn=zeros(length(f),1);	  % lorentz & birnbaum
lormix=zeros(length(f),1);  	  % lorentz & 1st order mixing
lormix_birn=zeros(length(f),1);   % lorentz & 1st order mixing & birnbaum

%%%%%%%%%%%%%%%%%% Compute lorentz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(freqr_shift)
  w_tot=(pressure_self*w_selfr(i)+pressure_for*w_forr(i))/pressure_ref;
  temp=strenrt(i)*(w_tot)*(f/freqr_shift(i))./((f-freqr_shift(i)).^2 +...
     (w_tot).^2);
  temp2=strenrt(i)*(f/freqr_shift(i)).*(w_tot+...
	   (f-freqr_shift(i))*ymix(i))./((f-freqr_shift(i)).^2+(w_tot).^2);

  lor=lor+temp;          % add up lorentz
  lormix=lormix+temp2;   % add up 1st order mixing

  chi=birnbaum(f',freqr_shift(i),w_tot,temperature,duration);  
  chi=chi';
  lor_birn=lor_birn+temp.*chi;          % add up lorentz & birnbaum
  lormix_birn=lormix_birn+temp2.*chi;   % add up 1st order mixing & birnbaum
  end

lor=lor*K_scale_lor*bsm;
lor_birn=lor_birn*K_scale_lor*bsm;
lormix=lormix*K_scale_lor*bsm;
lormix_birn=lormix_birn*K_scale_lor*bsm;

