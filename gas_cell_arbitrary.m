function [w,od] = gas_cell_arbitrary(iData,units,press,partpress,temperature,CellLength,iGas,f1,f2);

%% iData = +1 if data is coming in profile structure
%% iData = -1 if data is going to be interactively input
%% units = 1,2,3 == atm,torr,mb

%% based on gas_cell_water.m  gas_cell_co2.m  gas_cell_others.m

tc2k = 273.15;
MGC = 8.314674269981136;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 

if iData == -1
  spectral = input('Enter [start stop] wavenumber : ');

  f1 = spectral(1);
  f2 = spectral(2);

  units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 

  press       = input('Enter total pressure  : ');  
  temperature = input('Enter temperature (in C) : ');  

  partpress   = input('Enter gas partial pressure  : ');  
  CellLength  = input('Enter path cell length (in cm) : ');  

  iGas        = input('Enter gasID : ')
end

temperature = temperature + tc2k;

if iGas == 1
  addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
  [p_wat, p_ice]=svp(temperature);
  if units == 1
    %% change to atm
    p_wat =   p_wat/1013;
  elseif units == 2
    %% change to atm, then to torr
    p_wat =   p_wat/1013;
    p_wat = p_wat * 760;
  end
  fprintf(1,'at a temp of %8.4f K, the svp_wat is %8.6f \n',temperature+273,p_wat);
end

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
GasAmt = CellLength*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2   
gascellparams = [press partpress temperature GasAmt]; 

fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 

IPfname = ['IPFILES/gas_cell_arb_gasID_' num2str(iGas,'%02d')];
fprintf(1,'opening file %s \n',IPfname);

fid = fopen(IPfname,'w');
fprintf(fid,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams);
fclose(fid);

if iGas == 1
  disp(' runnning WATER ....')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% do run8water

  % divide   self    for    divide by                              result
  % ----------------------------------------------------------------------
  %   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF)
  %    1     1       0     q v tanh(c2 v/2T) (296/T) * ps           CS
  %    1     0       1     q v tanh(c2 v/2T) (296/T) * (p-ps)       CF
  %    1     ?       ?     q v tanh(c2 v/2T) (296/T)                ps CS + pf CF

  if f1 >= 500 & f1 < 2830
    toptsx = runXtopts_params_smart(2000);
    [w,d8kc]    = run8water(1,f1,f2,IPfname,toptsx);
    toptsxc.CKD  = 1
    toptsxc.ffin = toptsx.ffin;
    [wc,d8conkc] = run8watercontinuum(1,f1,f2,IPfname,toptsxc);
    d8conkc = boxint2(d8conkc,5);
  else
    [w,d8kc]    = run8water(1,f1,f2,IPfname);
    topts.CKD   = 1;
    [w,d8conkc] = run8watercontinuum(1,f1,f2,IPfname,topts);
    d8conkc = boxint2(d8conkc,5);
  end

  %disp('saving to gas_cell_water.mat ....');
  %str = input('Enter comment : ');
  %gasID = [1 101 102 103];
  %array = gascellparams;
  %saver = ['save gas_cell_water.mat str w d8kc d8conkc  gasID array']; eval(saver);

  whos d8conkc d8kc
  od  = d8conkc+d8kc;

  clf;
  plot(w,(d8conkc+d8kc))
  %fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 

elseif iGas == 2
  disp(' runnning SIMPLE CO2 ....')
  if f1 >= 500 & f1 < 2830
    toptsx = runXtopts_params_smart(2000);
    [w,od] = run8(2,f1,f2,IPfname,toptsx);
  else
    [w,od] = run8(2,f1,f2,IPfname);
  end

elseif iGas > 2
  disp(' runnning OTHER GAS ....')
  if f1 >= 500 & f1 < 2830
    toptsx = runXtopts_params_smart(2000);
    [w,od] = run8(iGas,f1,f2,IPfname,toptsx);
  else
    [w,od] = run8(iGas,f1,f2,IPfname);
  end
end