% kcartaoutput file from kCARTA has "atm" units
%units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 
units = 1;

tc2k = 273.15;
MGC = 8.314674269981136  ;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 

theprofile = load('CRIS/kcartawateroutput.txt');
%press       = input('Enter total pressure  : ');  
%partpress   = input('Enter gas partial pressure  : ');  
%temperature = input('Enter temperature (in K) : ');  
%GasAmt      = input('Enter gas amt (mol/cm2) ');  
%% have dumped out stuff for layer 19
whos theprofile
thegasid    = theprofile(:,2);
if mean(thegasid ~= 1)
  error('have not read in profiles for gasid = 1')
  end
theprofile  = theprofile(10,3:6);
press       = theprofile(1)
partpress   = theprofile(2)
temperature = theprofile(3)
GasAmt      = theprofile(4)

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
 
% no need to change to kmoles cm-2 as the warning.msg file from kCARTA has correct units
% GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2   

gascellparams = [press partpress temperature GasAmt]; 

fprintf(1,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams) 

fid = fopen('IPFILES/cris_cell_water_kcarta','w');
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
[w,d8kc]    = run8water(1,x1,x2,'IPFILES/cris_cell_water_kcarta');
topts.CKD = 1;
[w,d8conkc] = run8watercontinuum(1,x1,x2,...
                           'IPFILES/cris_cell_water_kcarta',topts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read kcarta product
[dkc,wkc] = readkcstd('CRIS/test_water10.dat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/sergio/SPECTRA
save CRIS/run8_kcartawater10  w d8kc d8conkc gascellparams dkc wkc

clf
subplot(211); plot(wkc,dkc(:,1),wkc,d8kc); grid; title('water lines')
subplot(212); plot(wkc,dkc(:,1)./d8kc');   grid

clf
subplot(211); plot(wkc,sum(dkc(:,[2:3])'),wkc,d8conkc); grid; 
  title('water continuum')
subplot(212); plot(wkc,sum(dkc(:,[2 3])')./d8conkc);   grid
