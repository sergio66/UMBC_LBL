% kcartaoutput file from kCARTA has "atm" units
%units = input('Enter pressure units (1) atm (2) torr (3) mb : '); 
units = 1;

tc2k = 273.15;
MGC = 8.314674269981136  ;    
torr2mb  = 1013.25 / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1/1013.25; 

theprofile = load('CRIS/kcartaoutput.txt');
%press       = input('Enter total pressure  : ');  
%partpress   = input('Enter gas partial pressure  : ');  
%temperature = input('Enter temperature (in K) : ');  
%GasAmt      = input('Enter gas amt (mol/cm2) ');  
%% have dumped out stuff for layer 57
whos theprofile
thegasid    = theprofile(:,2);
if mean(thegasid ~= 2)
  error('have not read in profiles for gasid = 2')
  end
theprofile  = theprofile(57,3:6);
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

fid = fopen('IPFILES/cris_cell_co2_kcarta','w');
fprintf(fid,'1  %12.8f  %12.8f   %6.3f   %10.5e \n',gascellparams);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% do run7co2 
[w,d7kc]    = run7co2(2,630,830,'IPFILES/cris_cell_co2_kcarta');
[w,d7simplekc] = run7(2,630,830,'IPFILES/cris_cell_co2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read kcarta product
[dkc,wkc] = readkcstd('CRIS/test_cris50.dat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% now also do GLAB run

seder = ['!sed -e "s/PPP/' num2str(press) '/g" -e "s/PXPXPX/' num2str(partpress) '/g"'];
seder = [seder ' -e "s/TTT/' num2str(temperature) '/g" -e "s/QQQ/' num2str(GasAmt) '/g" GLAB/Glab/glab_template.ip > GLAB/Glab/co2_kcarta.ip'];
eval(seder)
eval(['!cat GLAB/Glab/co2_kcarta.ip']);
cd GLAB/Glab
rmer = ['!/bin/rm co2_kcarta.rad  co2_kcarta.tau  co2_kcarta.txt  fort.30  fort.4  junk.dat'];
eval(rmer);
glaber = ['!run_avg co2_kcarta']; eval(glaber)
glab = load('GLAB/Glab/co2_kcarta.tau');

cd /home/sergio/SPECTRA
save CRIS/run7_kcarta50  w d7kc d7simplekc glab gascellparams dkc wkc