units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 

spectral = input('Enter [start stop] wavenumber : ');

tc2k = 273.15;
MGC = 8.314674269981136  ;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 

press       = input('Enter total pressure  : ');  
temperature = input('Enter temperature (in C) : ');  

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
[p_wat, p_ice]=svp(temperature+273);
if units == 1
  p_wat =   p_wat/1013;
elseif units == 2
  p_wat =   p_wat/1013;
  p_wat = p_wat * 760;
end
fprintf(1,'at a temp of %8.4f K, the svp_wat is %8.6f \n',temperature+273,p_wat);

partpress   = input('Enter gas partial pressure  : ');  
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
 
% change to kmoles cm-2
GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2   
gascellparams = [press partpress temperature GasAmt]; 

fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 

fid = fopen('IPFILES/gas_cell_water','w');
fprintf(fid,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% do run8water

% divide   self    for    divide by                              result
% ----------------------------------------------------------------------
%   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF)
%    1     1       0     q v tanh(c2 v/2T) (296/T) * ps           CS
%    1     0       1     q v tanh(c2 v/2T) (296/T) * (p-ps)       CF
%    1     ?       ?     q v tanh(c2 v/2T) (296/T)                ps CS + pf CF

x1 =  830; x2 =  930;
x1 = 2080; x2 = 2205;
x1 = spectral(1); x2 = spectral(2);

[w,d8kc]    = run8water(1,x1,x2,'IPFILES/cris_cell_water_kcarta');
topts.CKD = 1;
[w,d8conkc] = run8watercontinuum(1,x1,x2,...
                           'IPFILES/cris_cell_water_kcarta',topts);

disp('saving to gas_cell_water.mat ....');
str = input('Enter comment : ');
gasID = [1 101 102 103];
array = gascellparams;

saver = ['save gas_cell_water.mat str w d8kc d8conkc  gasID array']; eval(saver);

clf;
plot(w,(d8conkc+d8kc))
fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 
