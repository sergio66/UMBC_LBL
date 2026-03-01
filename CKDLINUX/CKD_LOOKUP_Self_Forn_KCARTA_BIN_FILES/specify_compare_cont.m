iA = input('Enter CKD vA (21,23,24,60 or 1,2,3,4,5,6 : ');
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' num2str(iA) '.bin'];
[ksA, freq, temp] = contread(fname);
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor' num2str(iA) '.bin'];
[kfA, freq, temp] = contread(fname);

iB = input('Enter CKD vB (21,23,24,60 or 1,2,3,4,5,6 : ');
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf' num2str(iB) '.bin'];
[ksB, freq, temp] = contread(fname);
fname = ['/home/sergio/KCARTADATA/General/CKDieee_le/CKDFor' num2str(iB) '.bin'];
[kfB, freq, temp] = contread(fname);

ii = find(temp == 300);

figure(1)
semilogy(freq,ksA(ii,:),freq,ksB(ii,:),'r','LineWidth',2);
legend(num2str(iA),num2str(iB)); xlabel('Wavenumber cm-1'); ylabel('CS(300K)')
axis([650 2800 5e-28 1e-23]); grid

figure(2)
semilogy(freq,kfA(ii,:),freq,kfB(ii,:),'r','LineWidth',2);
legend(num2str(iA),num2str(iB)); xlabel('Wavenumber cm-1'); ylabel('CF(300K)')
axis([650 2800 1e-31 1e-24]); grid;
