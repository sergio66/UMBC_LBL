% script file to produce BINARY fortran look up table at 1.0 cm-1 * X
% for bands other than thermal IR

disp('look at /home/sergio/SPECTRA/xckd_lookupBIN_IR.m');
disp('look at /home/sergio/SPECTRA/CKDLINUX/CKD_LOOKUP_ALLBANDS_May2013/ckd_lookupBIN_IR.m');
disp('ret to continue')
pause

disp('look at /home/sergio/SPECTRA/CKDLINUX/CKD_LOOKUP_ALLBANDS_Oct2014/ckd_lookupLOOP.m');
disp('ret to continue')
pause

cd /home/sergio/SPECTRA/CKDLINUX/CKD_LOOKUP_ALLBANDS_Oct2014/
foutROOT = '/asl/data/kcarta/H2012.ieee-le/';
main_ckd_lookupLOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
compare_all_bands_ckd_lookupBIN_ieee_le
%}

disp('can do eg ls -lt /asl/data/kcarta/H2012.ieee-le/*/etc.ieee-le/CKD*25* | more')
figure(1); clf
[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR15_30/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR15_30/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR30_50/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR30_50/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR50_80/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR50_80/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR80_150/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR80_150/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR140_310/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR140_310/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR300_510/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR300_510/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR500_605/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR500_605/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),freq,kf(21,:),'r'); hold on

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


