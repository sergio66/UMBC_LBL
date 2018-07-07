units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 
 
tc2k = 273.15;
MGC = 8.314674269981136  ;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 
press       = input('Enter total pressure  : ');  
partpress   = input('Enter gas partial pressure  : ');  
temperature = input('Enter temperature (in C) : ');  
GasAmt      = input('Enter path cell length (in cm) ');  
 
temperature = temperature + tc2k;

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
GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2   
array = [press partpress temperature GasAmt]; 

fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',array) 

fid = fopen('IPFILES/cris_cell_co2','w');
fprintf(fid,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',array);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d7]    = run7co2(2,630,830,'IPFILES/cris_cell_co2');
[wa,d8ha] = run8co2(2,630,730,'IPFILES/cris_cell_co2');   %% hartman
[wb,d8hb] = run8co2(2,730,830,'IPFILES/cris_cell_co2');   %% hartman
d8 = [d8ha d8hb]

[w,d7simple] = run7(2,630,830,'IPFILES/cris_cell_co2');

save run7_hartman w d7 d8