%%%%%%%%%%%%%%%%%%%%%%%%%%  set_fit.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file updates the values of the fitted parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bsm = flag4(1);
beta = b_fit(1);
duration = b_fit(2);
frequency_shift = b_fit(3);
c0 = b_fit(4);
cspan = b_fit(5);

band_strength_multiplier=bsm;

if pressure_for==0
      beta_pure=beta;
      beta_for=0;
else
      beta_for=((pressure_self+pressure_for)*beta-beta_pure*pressure_self)/...
		pressure_for;
end
