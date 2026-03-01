base = '/home/sergio/KCARTADATA/General/CKDieee_le/CKD';

fname = [base 'Self00.bin'];
[kxS00, freq, temp] = contread(fname);
fname = [base 'For00.bin'];
[kxF00, freq, temp] = contread(fname);

fname = [base 'Self21.bin'];
[kxS21, freq, temp] = contread(fname);
fname = [base 'For21.bin'];
[kxF21, freq, temp] = contread(fname);

fname = [base 'Self23.bin'];
[kxS23, freq, temp] = contread(fname);
fname = [base 'For23.bin'];
[kxF23, freq, temp] = contread(fname);

fname = [base 'Self24.bin'];
[kxS24, freq, temp] = contread(fname);
fname = [base 'For24.bin'];
[kxF24, freq, temp] = contread(fname);

fname = [base 'Self51.bin'];
[kxS51, freq, temp] = contread(fname);
fname = [base 'For51.bin'];
[kxF51, freq, temp] = contread(fname);

fname = [base 'Self55.bin'];
[kxS55, freq, temp] = contread(fname);
fname = [base 'For55.bin'];
[kxF55, freq, temp] = contread(fname);

fname = [base 'Self60.bin'];
[kxS60, freq, temp] = contread(fname);
fname = [base 'For60.bin'];
[kxF60, freq, temp] = contread(fname);

ii = 21;   %%temp(21) = 300
array = [kxS00(ii,:); kxS21(ii,:); kxS23(ii,:); kxS24(ii,:); ...
         kxS51(ii,:); kxS55(ii,:); kxS60(ii,:)];
semilogy(freq,array,'LineWidth',2);
grid; legend('00','21','23','24','51','55','60'); title('CS at 300 K')

ii = 17;   %%temp(17) = 260
array = [kxS00(ii,:); kxS21(ii,:); kxS23(ii,:); kxS24(ii,:); ...
         kxS51(ii,:); kxS55(ii,:); kxS60(ii,:)];
semilogy(freq,array,'LineWidth',2);
grid; legend('00','21','23','24','51','55','60'); title('CS at 260 K')