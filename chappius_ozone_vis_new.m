function ozoneOUT = chappius_ozone_vis(w,q);

% function ozone = chappius_ozone_vis(w,q);
% outputs an approzimation to O3 absorption between 400-800 nm
% based on (w) = input wavevector and (q) = input # of molecules/cm2
%%% this is based on the ozone_chappius.pdf paper
%%% see figure 3

%% gas cell length = 10 cm, p = 655.5 mmHg, T = 300K

MGC   = 8.314674269981136  ;  
press = 65.5/76.0;   %% p in atm
partpress = press; %% completely filled with ozone
temperature = 300;
GasAmt = 10;             %% input('Enter path cell length (in cm) ');
%change to kilomolcles cm-2 
GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmolec/cm2 

%% figure (3) has approx a triangle, from 4200 A to 7250 A
%% and peak abs coeff = 0.055 cm-1 at 6000 A

wavelength = 4200:1:7500;  ozone = zeros(size(wavelength));
p1 = [4250 0];  pA = [6000 0.055]; p2 = [7500 0];
slope1 = 0.055/(6000-4250);        slope2 = 0.055/(7500-6000);

i1 = find(wavelength >= 4250 & wavelength <= 6000);
ozone(i1) = slope1*(wavelength(i1)-4250);

i2 = find(wavelength > 6000 & wavelength <= 7500);
ozone(i2) = slope2*(7500 - wavelength(i2));

plot(wavelength,ozone)

wavelength = wavelength/10000;
plot(wavelength,ozone); title('1'); xlabel('\lambda um'); pause(0.1)

wavenumber = 10000./wavelength;
plot(wavenumber,ozone); title('2'); xlabel('\nu cm-1'); pause(0.1)

%% the units of absorption here are cm-1
%% can change it to units of per molecule/cm2 by 
%%   multiplying by gas cell length  (get total absorption of all molecules)
%%   dividing by GasAmt

GasAmt = GasAmt*1000;   %% molecules/cm2
ozone = ozone*10/GasAmt;
plot(wavelength,ozone);  title('3');  xlabel('\lambda um'); pause(0.1)

%%% of course, then we rescale to 
%% GEOPHYSICAL RESEARCH LETTERS, VOL. 19, NO. 9, PAGES 933–936, 1992

%% Laser measurements of ozone absorption cross sections in the Chappuis band
%% Stuart M. Anderson
%% Department of Physics, Augsburg College
%% Konrad Mauersberger
%% School of Physics and Astronomy, University of Minnesota
%% Abstract

%% We have developed a sensitive spectrometer which exploits several He-Ne 
%% laser transitions in the visible for precise, high resolution measurements 
%% of Chappuis band ozone absorption cross sections at room temperature. An 
%% overall uncertainty of better than 1% has been achieved through a 
%% combination of transducer calibrations and an experimental procedure which 
%% unambiguously determines the impurity content ot each ozone sample. The 
%% cross section near 602.5 nm is (5.16±0.03)×10−21 cm2/molecule. Our results
%% compare favorably with those from most previous studies in the visible 
%% range, although some are clearly lower than this consensus probably due to 
%% errors in ozone density measurements. © American Geophysical Union 1992 

ozone = ozone/max(ozone) * 5.16e-21;
plot(wavelength,ozone);   xlabel('\lambda cm'); pause(0.1)

ozone = ozone * q;   %% find optical depth

ozoneOUT = zeros(size(w));
iX = find(w >= min(wavenumber) & w <= max(wavenumber));
if length(iX) > 0
  ozoneOUT(iX) = interp1(wavenumber,ozone,w(iX),[],'extrap');
  end
figure(1); plot(w(iX),ozoneOUT(iX)); xlabel('\nu cm-1'); ylabel('abs coeff')

lala= load('/home/sergio/SPECTRA/VISIBLE_OD/O3.csv');
T = temperature;

[min(lala(:,1)) max(lala(:,1))]

lala(:,1) = lala(:,1)/10000;   %%% wavelength

od = lala(:,2) + lala(:,3)*T + lala(:,4)*T.*T;

toto = interp1(10000./lala(:,1),od,w);
figure(2); semilogy(w,toto);

figure(1)