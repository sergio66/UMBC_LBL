% test parameters
gid = 2;       % HITRAN gas id 
pT  = 40;      % cell pressure in torr
tC  = 17;      % cell temperature in degrees C
Lcm = 12.59;   % path length in cm
v1  = 605;     % band start
v2  = 1230;    % band end  

v1  = 655;     % band start
v2  = 705;     % band end  

% convert pressures
torr2mb = 1013.25 / 760;
mb2atm = 9.86923267e-4;
patm = pT * torr2mb * mb2atm;
ppatm = patm;

% convert temperature
tC2K = 273.15;
tK = tC + tC2K;

% change to kmoles/cm2
MGC = 8.314674269981136;
gAmt = Lcm * 101325 * ppatm / 1e9 / MGC / tK;

% write a params file
gf = './run8_CO2.tmp';
fi = fopen(gf, 'w');
fprintf(fi, '%4d%12.4e%12.4e%8.2f%12.4e\n', 1, patm, ppatm, tK, gAmt);
fclose(fi);

% call run8
addpath /home/sergio/SPECTRA
[fr, absc] = run8(gid, v1, v2, gf);
[fr, abscL] = run8co2_linemixUMBC(gid, v1, v2, gf);
[fr, abscV] = run8co2_voigt(gid, v1, v2, gf);

% keep data in column order
fr = fr(:); absc = absc(:);

% save the results
save run8_40t_CO2 fr absc*
