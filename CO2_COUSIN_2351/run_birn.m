%%%%% where you see %%%%% means i have modified this slightly than from
%%%%% the  original version found in /beet/users/tobin/Co2pr/Fits_new/Fit_all

loader_mod 

%%%%eval(['load ../Back/' fil_nam '_back'])

bsm=input('Enter bsm: ');
band_strength_multiplier=bsm;
duration=input('Enter tau2: ');
%%%%beta_for=input('Enter beta(for): ');
%%%%beta_pure=input('Enter beta(CO2): ');
%%%%beta=(pressure_self*beta_pure+pressure_for*beta_for)./...
%%%%     (pressure_for+pressure_self);

beta=input('Enter beta: ');
beta_pure=input('Enter beta(CO2): ');

frequency_shift=input('Enter frequency_shift: ');
%%%%K_back=mix_birn;

%%%%keyboard

efitter
trans_pop
wfunco2er
y1ser
klormix_birn
%%%%full_mix4_2
%%%%%eval(['save ' fil_nam '_calc'])



