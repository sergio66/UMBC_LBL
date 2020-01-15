%% see run8.m

  MGC = 8.314674269981136  ;
  press       = input('Enter total pressure (in atm) : ');
  partpress   = input('Enter gas partial pressure (in atm) : ');
  temperature = input('Enter temperature (in K) : ');
  GasAmt      = input('Enter path cell length (in cm) ');
  %change to kmoles cm-2
  GasAmt = GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2

fprintf(1,'GasAmt = %10.6e  kmoles/cm2 \n',GasAmt)
