oo = load('output');  oo2 = load('output2'); 
dd = load('LM_TEST.dat');
whos
plot(oo(:,1)-dd(1:4500,1))
plot(oo(2:4500,1)-dd(1:4499,1))
plot(oo(3:4500,1)-dd(1:4498,1))
oo = oo(3:4500,:); oo2 = oo2(1:4498,:);  dd = dd(1:4498,:);
ii=1; plot(oo(:,1),oo(:,ii)./dd(:,ii))
ii=2:4; plot(oo(:,1),oo(:,ii)./dd(:,ii))
ii=2:4; plot(oo(:,1),oo(:,ii)./(dd(:,ii)*3.315e7))
ii=2:4; plot(oo(:,1),oo2(:,ii)./(dd(:,ii)))

ii=2:4; plot(oo(:,1),oo(:,ii)./(dd(:,ii)*3.315e7), oo(:,1),oo2(:,ii)./(dd(:,ii)))
ii=input('enter (1) L (2) ist order (3) Full : '); 
ii = ii + 1;
  plot(oo(:,1),oo(:,ii)./(dd(:,ii)*3.315e7), oo(:,1),oo2(:,ii)./(dd(:,ii))) 
  disp('ret to continue : '); pause
  ll = input('enter artificial length : ');   %% artificial length
  ll = abs(ll); ll = max(ll,1);
  plot(oo(:,1),exp(-ll*oo(:,ii))./exp(-ll*dd(:,ii)*3.315e7), oo(:,1),exp(-ll*oo2(:,ii))./exp(-ll*dd(:,ii)))

