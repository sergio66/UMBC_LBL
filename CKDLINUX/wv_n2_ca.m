function [ca] = wv_n2_ca(freq,Ptot,PN2,PWV,T)

% ref : Effect of humiity on the absorption continua of CO2 abd N2 near 4 um :
% calculations, comparisons with measurements, consequences on atmospheric spectra
% JM Hartmann, C. Boulet, D. Tran, H. Tran, Y. Baranov
% Journal of Chemical Phycis, v 148 pg 54304 (2018)
%
% this is a wee bit hard to do, since I have no idea how to do the N2/WV->N2/WV transitions
% just look at N2_routines.f

% Section III B Equation (2)
addpath /home/sergio/SPECTRA

%T = 300;

eta = [12.1 12.1 18.0 26.4];

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;

%freq = 2000:0.0025:2800;
%freq = 2000:0.1:2800;
ca = zeros(size(freq));

[iYes,line1] = findlines_plot(2000,3000,1,2016);
[iYes,line22] = findlines_plot(2000,3000,22,2016);

line = line22;
band_center = 2329;

for ii = 1 : line.linct
  deltasigma = (freq - band_center);
  deltasigma = (freq - line.wnum(ii));  
  G = zeros(size(freq));
  for kk = 1 : 1
    G = G + 2/3/pi/eta(kk)*(deltasigma/eta(kk)).^2./(1+exp(-c2/T*deltasigma)).*besselk(2,abs(deltasigma/eta(kk)));
  end
  junk = line.stren(ii) * G; %% need to adjust line stren for T
  ca = ca + junk;
end

alpha = 0.007297352566417;  %% https://en.wikipedia.org/wiki/Fine-structure_constant
a0    = 5.2917721067e-9;    %% in cm, https://en.wikipedia.org/wiki/Bohr_radius
no    = 2.6867805e+19;      %% amagat molecules/cm3

ca = ca*4*pi*pi/3*alpha*(a0^5)*no*no;
ca = ca*PN2/Ptot*PWV/Ptot;

plot(freq,ca)
