units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 

if ~exist('gasID')
  gasID = input('Enter gasID ');
end
spectral = input('Enter [start stop] wavenumber : ');
 
tc2k = 273.15;
MGC = 8.314674269981136  ;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 

press       = input('Enter total pressure  : ');  
partpress   = input('Enter gas partial pressure  : ');  
temperature = input('Enter temperature (in C) : ');  
GasAmt      = input('Enter path cell length (in cm) ');  
% str         = input('Enter comment : ');

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
GasAmt = GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2   
gascellparams = [press partpress temperature GasAmt]; 

fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 


if ~exist('tempfilename')
  tempfilename = '/home/sergio/git/SPECTRA/IPFILES/gas_cell';
  tempfilename = 'IPFILES/gas_cell';
end

profname = tempfilename;  
fid = fopen(tempfilename,'w');
fprintf(fid,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wx,dx]   = run8(gasID,spectral(1),spectral(2),tempfilename);

iSave = input('Save to gas_cell_co2.mat (-1 NO default/+1 YES)  : ');
if length(iSave) == 0
  iSave = -1;
end
if iSave > 0
  disp('saving to gas_cell.mat ....');
  save gas_cell.mat wx dx gasID gascellparams
end
w = wx;
d = dx;
