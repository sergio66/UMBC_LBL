dir0 = '/home/sergio/KCARTADATA/General/CKDieee_le/';
dirN = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/'

fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf24.bin';
[ks24, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor24.bin';
[kf24, freq, temp] = contread(fname);

%fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf60.bin';
%[ks60, freq, temp] = contread(fname);
%fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor60.bin';
%[kf60, freq, temp] = contread(fname);

fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf1.bin';
[ks1, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor1.bin';
[kf1, freq, temp] = contread(fname);

fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf4.bin';
[ks4, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor4.bin';
[kf4, freq, temp] = contread(fname);

ii = find(temp == 300);

figure(1)
semilogy(freq,[ks24(ii,:); ks1(ii,:); ks4(ii,:)],'LineWidth',2);
legend('24','1','4'); xlabel('Wavenumber cm-1'); ylabel('CS(300K)')
axis([650 2600 5e-28 1e-23]); grid

figure(2)
semilogy(freq,[kf24(ii,:); kf1(ii,:); kf4(ii,:)],'LineWidth',2);
legend('24','1','4'); xlabel('Wavenumber cm-1'); ylabel('CF(300K)')
axis([650 2600 1e-31 1e-24]); grid;

fid = fopen('ckd_compare_24_1_4.txt','w');
array = ...
  [freq; ks24(ii,:); ks1(ii,:); ks4(ii,:); kf24(ii,:); kf1(ii,:); kf4(ii,:)];
fprintf(fid,'%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e \n',array);
fclose(fid);

clear all
dd = load('ckd_compare_24_1_4.txt');
figure(3); semilogy(dd(:,1),dd(:,2:7))
axis([600 2800 1e-32 1e-20]); grid