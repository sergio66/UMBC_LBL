%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y1s.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This program computes the first order mixing coefficients: Y_1st.
%   
%   They are calculated with the equation:
%
%	Yi= 2* SUM{(dk/di)*Wki/(vi-vk)} 
%		where the sum is over k not equal to i
%		with dk,i= dipole matrix elements for lines k,i
%		     Wki = W matrix element for relaxation from i to k
%		     vk,i= wavenumber of line k,i
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y_1st=y1s(j,freq,elower,stren,W_plus,prb)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor;
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back
global K_back strenrt w_forr w_selfr

if prb=='R'
	trans_ampl=trans_amplr;
end
if prb=='P'
	trans_ampl=trans_amplp;
end

no_lines=length(freq);Y_1stmix=zeros(1,no_lines);
freq_shift=freq+frequency_shift/100;

W_plus2=W_plus.*(beta+(1-beta)*eye(no_lines,no_lines)); 
for n=1:no_lines
  for m=1:no_lines
    if m~=n
      Y_1stmix(n)=Y_1stmix(n)+trans_ampl(m)*(W_plus2(m,n))/(freq(n)-freq(m));
    end
  end
end
Y_1st=2*Y_1stmix./trans_ampl';
