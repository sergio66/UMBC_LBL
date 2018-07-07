path(path,'/home/sergio/SPECTRA');

global betafudge docfudge
global beta beta_pure beta_for bsm duration frequency_shift fudge
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor pick_no
global trans_amplr populationr population_tr t_rawdata;
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back K_back strenrt w_forr w_selfr abscoef 
global sigsig_comp sigsigR sigsigP ymixR wwwR 

choice1=2;
if choice1==1
   beta_for        = input('Enter beta(N2) : ');
   beta_pure       = input('Enter beta(CO2): ');
   beta            = input('Enter beta : ');
   bsm             = input('Enter bsm: ');
   duration        = input('Enter duration: ');
   frequency_shift = input('Enter freq_shift*100: ');
else
   beta_for=1;beta_pure=0.75;beta=0.95;
   bsm=1; band_strength_multiplier=bsm;
   duration=0.005; 
   frequency_shift=0; 
end

band_strength_multiplier=bsm;
loader_mod_COUSIN;
