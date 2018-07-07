%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  tdif_birn.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the difference between 1st order mix transmission
%       with birnbaum and actual transmission.
%
%  same as tdif_birn_Oct99 except that it gets rid of the magic numbers
%  in the Cousin-Mixing blending
%
% also it is modified so that if bsm < 0, or if duartion,beta,freq shift are 
% unacceptable, it gives 999999999999999999999999
% else it will go and happily call tdif_birn_goodparams.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function final_diff=tdif_birn(bstart,freqr,freqp,jp,f,flag4)

global beta beta_pure beta_for bsm duration frequency_shift fudge
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back pick_no
global K_back strenrt w_forr w_selfr abscoef 
global sigsig_comp sigsigR sigsigP ymixR wwwR

iOK = 1;     %%%assume things are ok

beta = bstart(1);
if ((beta < 0) | (beta > 1))
  iOK = -1;
  end

duration = bstart(2);
if ((duration < 1e-3) | (duration > 28e-3))
  iOK = -1;
  end

frequency_shift = bstart(3);
if (abs(frequency_shift) > 0.5)
  iOK = -1;
  end

c0 = bstart(4);
if c0 < 1.0
  iOK = -1;
  end

cspan = bstart(5);
if cspan < 0
  iOK = -1;
  end

band_strength_multiplier=flag4;

if (iOK > 0)
  final_diff=tdif_birn_goodparams(bstart,freqr,freqp,jp,f,flag4);
else
  final_diff=ones(size(f))*1e4;
  end