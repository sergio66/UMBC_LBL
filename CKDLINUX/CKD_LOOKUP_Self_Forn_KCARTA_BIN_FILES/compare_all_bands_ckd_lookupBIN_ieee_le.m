%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('can do eg ls -lt /asl/data/kcarta/H2012.ieee-le/*/etc.ieee-le/CKD*25* | more')
figure(1); clf
[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR15_30/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR15_30/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR30_50/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR30_50/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR50_80/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR50_80/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR80_150/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR80_150/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR140_310/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR140_310/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR300_510/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR300_510/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR500_605/etc.ieee-le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2012.ieee-le/FIR500_605/etc.ieee-le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor25.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf25.bin');
figure(1); semilogy(freq,ks(21,:),'b',freq,kf(21,:),'r','linewidth',4); hold on

hold off
title('blue = Self red = Forn')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
disp('can do eg ls -lt /asl/data/kcarta/*/*/etc.ieee-le/CKD* | more');
%[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR15_30/etc.ieee-le/CKDFor1.bin');
%[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR15_30/etc.ieee-le/CKDSelf1.bin');
%figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

%[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR30_50/etc.ieee-le/CKDFor1.bin');
%[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR30_50/etc.ieee-le/CKDSelf1.bin');
%figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

%[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR50_80/etc.ieee-le/CKDFor1.bin');
%[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR50_80/etc.ieee-le/CKDSelf1.bin');
%figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR80_150/etc.ieee-le/CKDFor1.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR80_150/etc.ieee-le/CKDSelf1.bin');
figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR140_310/etc.ieee-le/CKDFor1.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR140_310/etc.ieee-le/CKDSelf1.bin');
figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR300_510/etc.ieee-le/CKDFor1.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR300_510/etc.ieee-le/CKDSelf1.bin');
figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR500_605/etc.ieee-le/CKDFor1.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/H2004_otherbands.ieee-le/FIR500_605/etc.ieee-le/CKDSelf1.bin');
figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

[kf,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor1.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf1.bin');
figure(1); semilogy(freq,ks(21,:),'c',freq,kf(21,:),'m','linewidth',2); hold on

disp('CKD 24')
[kf,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor24.bin');
[ks,freq,T] = contread2('/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf24.bin');
figure(1); semilogy(freq,ks(21,:),'g',freq,kf(21,:),'k','linewidth',2); hold on

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on



