function [CA] = wv_co2_ca(T,p,v);

%% GEISA to HITRAN http://eodg.atm.ox.ac.uk/RFM/geihit.html
%% GEISA at http://cds-espri.ipsl.upmc.fr/etherTypo/?id=1732

%% T = temperature (K)
%% p = pressure (atm)
%% v = vector of wavenumbers = vector

%% now sum these over individual lines, pretend this is the 2400 cm-1 bandhead line

%load ../CO2_MATFILES/hit2350.mat
lines = load('../CO2_MATFILES/hit2350.mat');

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;
c22 = c2/2;

c2T  = c2/T;
c22T = c22/T;

%see qtips16.m
    abcd = [ ... 
  -1.65407e+00 , 9.58162e-01, -7.65288e-04 , 2.74493e-06 ; 
  -2.44629e+00 , 1.90180e+00, -1.49188e-03 , 5.66308e-06 ; 
  -3.55661e+00 , 2.03289e+00, -1.65085e-03 , 5.94879e-06 ; 
  -2.06391e+01 , 1.18565e+01, -9.55268e-03 , 3.43467e-05 ; 
  -5.24975e+00 , 4.03437e+00, -3.21680e-03 , 1.22759e-05 ; 
  -3.04229e+01 , 2.35296e+01, -1.86113e-02 , 7.08549e-05 ; 
  -1.91991e+00 , 1.08022e+00, -8.92777e-04 , 3.23212e-06 ; 
  -2.21912e+01 , 1.25883e+01, -1.03145e-02 , 3.72636e-05 ; 
  -6.43127e+01 , 3.66937e+01, -2.98289e-02 , 1.07509e-04 ; 
  -2.82091e+00 , 2.14310e+00, -1.73806e-03 , 6.66988e-06 ; 
  -3.25966e+01 , 2.49767e+01, -2.00829e-02 , 7.68876e-05 ]; 
%  disp('The isotopes are 626  636  628  627  638  637  828  827  727  838  837'); 
a = abcd(1,1);
b = abcd(1,2);
c = abcd(1,3);
d = abcd(1,4);
qfcn = q(a,b,c,d,0,lines,T);
strengthT = find_stren(qfcn',v,T,lines.elower,lines.stren,1.00);
strengthT = strengthT/6.022045e26; %% find_stren multiplies line strength by kAvog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T0 = 273.15;
rbranch = find(lines.j_lower(:,5) == 'R');
rbranch = 1:length(lines.freq);
CA = zeros(size(v));
for iix = 1 : length(rbranch)
  ii = rbranch(iix);
  junk = zeros(size(v));
  
  blah1 = strengthT(ii) * exp(c22T*(v-lines.freq(ii))) .* (1-exp(-c2T*v))./(1-exp(-c2T*lines.freq(ii)));

  gammaW_CO2 = lines.w(ii)/10*p;
  gammaW_CO2 = lines.w(ii)*T/T0;

  % this is from the Can J. Phys paper
  % gammaW_CO2 = T/T0 * ones(size(lines.w(ii)))*0.13;
  
  chi = wv_co2_chi(T,abs(v-lines.freq(ii)));  
  blah2 = gammaW_CO2 * chi./abs(v-lines.freq(ii))./abs(v-lines.freq(ii));
  
  junk = blah1 .*(v/lines.freq(ii)) .* blah2;
  
  bad = find(abs(lines.freq(ii)-v) < 5);
  junk(bad) = 0;
  CA = CA + junk/pi;
end

CA = CA * 2.68678e19;     %% HITRAN strength = cm-1/(molecules/cm2) --- multiply by 1 amagat to convert to cm-2/amagat-1
semilogy(v,CA); grid
xlabel('Wavenumber cm-1'); ylabel('CA_{CO2-H2O} (cm-1/amagat^2)')
set(gca,'fontsize',16)

%{

values of -3.62(14)×10 -3 and -3.33(16)×10 -3 cm -1 /atm at 294 K are
obtained for O (2) and O (3), respectively.  For an easier comparison
with li terature values, we have converted these values in cm -1
/amagat which is the unit most frequently used for H 2 line shift (1
amagat is the molecular density corresponding to an ideal gas at 1 atm
and 0°C).  It leads to the values - 3.90(15)×10 -3 and -3.59(17)×10 -3
cm -1 /amagat for O (2) and O (3), respectively. T
%}
