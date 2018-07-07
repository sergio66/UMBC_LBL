%%%%%%%%%%%%%%%%%%%%%%%%%%  tdifer_birn.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file simply calls tdif_birn.m to fit for beta, band strength 
%	multiplier, constant frequency shift.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back pick_no
global K_back 
global w_forr w_selfr jp

options(1)=1;

for i=1:length(temperature)
  eval(['global W_plus_r' num2str(i)])
  eval(['global W_plus_p' num2str(i)])
  end

b_fit=leastsq('tdif_birn',bstart,options,[],freqr,freqp,jp,f,flag4);

