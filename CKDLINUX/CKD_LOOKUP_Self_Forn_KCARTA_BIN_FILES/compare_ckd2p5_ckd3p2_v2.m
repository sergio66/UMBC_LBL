addpath /home/sergio/SPECTRA
dirN = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('reading in raw files supplied with CKD')
ckdraw25 = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD2.5/cntnm/WATER.COEF_matlab_296K');
ckdraw25 = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD2.5/cntnm/WATER.COEF_matlab_300K');

ckdraw32 = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/run_example/WATER.COEF_matlab_296K');
ckdraw32 = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/run_example/WATER.COEF_matlab_300K');

figure(3); plot(ckdraw25(:,1)./ckdraw32(:,1))

figure(3); ii=2:3; plot(ckdraw25(:,1),ckdraw32(:,ii)./ckdraw25(:,ii))
  title('raw CKD32/CKD25'); hl = legend('self','forn');
figure(4); ii=2:3; plot(ckdraw25(:,1),ckdraw32(:,ii)./ckdraw25(:,ii))
  axis([600 3000 0.5 2])
  title('raw CKD32/CKD25'); hl = legend('self','forn');
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
%wax = load('xWATER.COEF');
wax = ckdraw25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CKD 1 or 6
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf6.bin';
[ks6, freq6, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor6.bin';
[kf6, freq6, temp] = contread(fname);
t300 = find(temp >= 300,1);

semilogy(wax(:,1),wax(:,2),'b',wax(:,1),wax(:,3),'r',...
         freq6,ks6(t300,:),'c',freq6,kf6(t300,:),'m')
semilogy(wax(:,1),wax(:,2),'b',wax(:,1),wax(:,3),'r',...
         freq6,ks6(1,:),'c',freq6,kf6(1,:),'m')
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CKD 25
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf25.bin';
[ks25, freq25, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor25.bin';
[kf25, freq25, temp] = contread(fname);
t300 = find(temp >= 300,1);

semilogy(wax(:,1),wax(:,2),'b',wax(:,1),wax(:,3),'r',...
         freq25,ks25(t300,:),'c',freq25,kf25(t300,:),'m')
semilogy(wax(:,1),wax(:,2),'b',wax(:,1),wax(:,3),'r',...
         freq25,ks25(1,:),'c',freq25,kf25(1,:),'m')
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
topts.useruser = +1;
topts.CKD      = 25;
[fr,krun_25] = run8watercontinuum(1,500,3000,'IPFILES/junk',topts);
topts.CKD      = 32;
[fr,krun_32] = run8watercontinuum(1,500,3000,'IPFILES/junk',topts);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CKD 32
%% first field is wavenumber cm-1
%% next  two are Cs and Cf without radiation field  (used by kCARTA)     1/(cm-1 molec/cm**2)
%% final two are Cs and Cf with    radiation field  (not used by kCARTA) 1/(molec/cm**2)
%%
woo = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/run_example/woo296');
woo = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/run_example/woo300');

fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf32.bin_orig';
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf32.bin';
[ks32, freq32, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor32.bin_orig';
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor32.bin';
[kf32, freq32, temp] = contread(fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t300 = find(temp >= 300,1);
figure(1); semilogy(freq25,ks25(t300,:),'b',freq32,ks32(t300,:),'r',freq6,ks6(t300,:),'k'); 
  title('Self (b)2.5 (r) 3.2 (k)1,6')
figure(2); semilogy(freq25,kf25(t300,:),'b',freq32,kf32(t300,:),'r',freq6,kf6(t300,:),'k'); 
  title('Forn (b)2.5 (r) 3.2 (k)1,6')
disp('ret to continue'); pause

%figure(1); semilogy(freq,ks25,'b',freq,kf32,'r'); title('Self (b)2.5 (r) 3.2')
%figure(2); semilogy(freq,kf25,'b',freq,ks32,'r'); title('Forn (b)2.5 (r) 3.2')
%figure(2); semilogy(freq,ks32,'b',freq,kf32,'r'); title('Forn (b)2.5 (r) 3.2')

plot(freq25,ks25,'b',freq32,ks32,'r'); title('Self')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f10 = 600:10:3000;
for ii = 1 : 31
  kss25(ii,:) = interp1(freq25,ks25(ii,:),f10);
  kff25(ii,:) = interp1(freq25,kf25(ii,:),f10);
  kss32(ii,:) = interp1(freq32,ks32(ii,:),f10);
  kff32(ii,:) = interp1(freq32,kf32(ii,:),f10);
end

plot(f10,kss32./kss25,'b',f10,kff32./kff25,'r')
title('from kCARTA (b)s32/s25 (r) f32/f25')

ksraw25 = interp1(ckdraw25(:,1),ckdraw25(:,2),f10);
kfraw25 = interp1(ckdraw25(:,1),ckdraw25(:,3),f10);
ksraw32 = interp1(ckdraw32(:,1),ckdraw32(:,2),f10);
kfraw32 = interp1(ckdraw32(:,1),ckdraw32(:,3),f10);

figure(1)
plot(f10,kss32(t300,:)./kss25(t300,:),'b',f10,kff32(t300,:)./kff25(t300,:),'r',...
     f10,ksraw25./kss25(t300,:),'c.-',f10,kfraw25./kff25(t300,:),'m.-');
title('roughly 296-300 K from kCARTA and raw')
  hl = legend('s32/s25','f32/f25','raws25/s25','rawf25/f25'); grid

figure(2)
plot(f10,kss32(t300,:)./kss25(t300,:),'b',f10,kff32(t300,:)./kff25(t300,:),'r',...
     f10,ksraw32./kss32(t300,:),'c.-',f10,kfraw32./kff32(t300,:),'m.-');
title('from kCARTA and raw')
  hl = legend('s32/s25','f32/f25','raws32/s32','rawf32/f32'); grid

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% with no radiation term ie just like kCARTA this needs v tanh*v/T)
figure(1)
semilogy(woo(:,1),woo(:,2),'b',freq25,ks25(t300,:),'r')
axis([500 3000 1e-27 1e-23]); grid on; title('cs')
hl = legend('2.5','3.2','location','best'); set(hl,'fontsize',10);

figure(2)
semilogy(woo(:,1),woo(:,3),'b',freq25,kf25(t300,:),'r')
axis([500 3000 1e-31 5e-25]); grid on; title('cf')
hl = legend('2.5','3.2','location','best'); set(hl,'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% from od300, these are molecular amounts  (MOL/CM**2)
%%   H2O       CO2        O3       N2O        CO       CH4        O2      OTHER
%% 2.421E+17 8.351E+15 0.000E+00 0.000E+00 0.000E+00 0.000E+00 5.083E+18 1.910E+19

P = 1013.0;  %% mb
T = 300.0;   %% K
L = 1.0;     %% cm
WVMR = 0.01; %% WV vmr

rho = (P*100)/8.3144598/T;   %% moles/m3
rho = rho*6.022140857e23;    %% molecules/m3
rho = rho/1e6;               %% molecules/cm3

Qtotal = rho*L;              %% molecules/cm2

QWV = WVMR*Qtotal*(1-WVMR);  %% have to do it wrt dryair hence Qtotal*(1-WVMR) = amount of dryair
Qdry = Qtotal*(1-WVMR);

%%%%%%%%%%%%%%%%%%%%%%%%%

% from MT_CKD3.2/cntnm/src/cntnm_progr.f

PAVE = P; TAVE = T; xlength = L; VMRH2O = WVMR;

XLOSMT = 2.68675E+19;   %% number of molecules at at 273 K, 1 atm
WTOT = XLOSMT*(PAVE/1013.)*(273./TAVE)* xlength;
W_dry = WTOT * (1.-VMRH2O);
WA     = 0.009     * W_dry;  %% argon
WN2    = 0.78      * W_dry;  %% N2
WK(7)  = 0.21      * W_dry;  %% O2
WK(2)  = 345.E-06  * W_dry;  %% CO2

if (abs(VMRH2O-1.) < 1.e-05) then
  wk(1) = wtot;
else
  WK(1) = VMRH2O * W_dry;
end

WBROAD=WN2+WA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% self continumm scales as Rself/Ro
%% forn continuum scales as (Rtot-Rself)/Ro
%% where R is density R = P/Po * To/T

v = woo(:,1);
% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;

radfield = v .* tanh(c2/2*v/T);
plot(v,woo(:,2).*radfield ./ woo(:,4))  %% should be 1
plot(v,woo(:,3).*radfield ./ woo(:,5))  %% should be 1

cself = QWV*P/1013*WVMR*woo(:,2).*radfield;
cforn = QWV*P/1013*(1-WVMR)*woo(:,3).*radfield;
odsergio = cself + cforn;

%{
odlblrtm = load('/home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/run_example/od300');
figure(3); plot(odlblrtm(:,1),odlblrtm(:,2),'b.-',v,odsergio,'r',fr,krun_25/15,'ko-',,fr,krun_32/15,'go-')
axis([fr(1) fr(end) 0 2e-4]); grid
%}
