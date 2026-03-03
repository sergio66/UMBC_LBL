units = input('Enter pressure units (1) atm (2) torr (3) mb : ');

%% run7* needs ATM, ATM, K, kmoles cm-2
set_c1_c2_avog_pref_tref

torr2mb     = (PREF/100)  / 760;
torr2atm    = 1 / 760;
mb2atm      = 1/(PREF/100);
press       = input('Enter total pressure (default atm) : '); 
partpress   = input('Enter gas partial pressure (default atm) : '); 
temperature = input('Enter temperature (in K) : '); 
GasAmt      = input('Enter path cell length (in cm) '); 

if units == 1
  press  = press;
  partpress = partpress;
elseif units == 2
  press  = press * torr2atm;
  partpress = partpress * torr2atm;
elseif units == 3
  press  = press * mb2atm;
  partpress = partpress * mb2atm;
  end

%change to kmoles cm-2  
GasAmt = GasAmt * PREF * partpress/1e9/MGC/temperature; %change to kmoles/cm2  
array  = [press partpress temperature GasAmt];
fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',array)
