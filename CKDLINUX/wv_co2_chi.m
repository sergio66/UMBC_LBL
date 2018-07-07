%% see Icarus v306 (2018) og 116-121 : Tran, Turbet, Chelin, Landsheere
%% MEasurements and modeling of absorption by CO2+H2O cmixtures in the spectral region beyond the CO2 v3 band

%% see KCARTA/SRCv1.18/kcont_xsec.f  ...
%% subr f_by_T_CKDfilereader --> ComputeCKD_Linear_March2002 --> aux_ckd

function [chi] = wv_co2_chi(T,dv0)

%% T = temperature = scalar
%% dv = vi - v where vi = linecenter, v = vector of wavenumbers = vector

v1 = 5;
v2 = 35;
v3 = 170;

B1 = 0.0689 - 2.4486./T + 64.085./(T.^2);
B2 = 0.00624 + 3.7273 ./ T - 299.144./T./T;
B3 = 0.0025;

dv = abs(dv0);

ii = find(dv <= v1);           chi(ii) = 1.0;
ii = find(dv > v1 & dv <= v2); chi(ii) = exp(-B1*(dv(ii)-v1));
ii = find(dv > v2 & dv <= v3); chi(ii) = exp(-B1*(v2-v1) - B2*(dv(ii)-v2));
ii = find(dv > v3);            chi(ii) = exp(-B1*(v2-v1) - B2*(v3-v2) - B3*(dv(ii)-v3));

%figure(1); plot(dv+2400,chi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
