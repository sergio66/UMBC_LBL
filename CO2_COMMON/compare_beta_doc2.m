ido = input('do the RAL (+1) or JJOHNS (-1) p/ps parameters : ');

%%%region == 1 ==> 2380-2391
%%%region == 2 ==> 2391-2400
%%%region == 3 ==> 2400-2405

if ido > 0
  %%%these are the RAL params
  pt=[2.6800e+01   1.6140e+02   5.6030e+02   9.6170e+02];  pt=pt/1013.25;
  ps=[2.6800e+01   2.6800e+01   2.6800e+01   2.6800e+01];  ps=ps/1013.25;
else
  %%%these are the JJOHNS params
  pt=[2.6800e+01 1.0153e+03 1.0033e+03];  pt=pt/1013.25; 
  ps=[2.6800e+01 5.0716e+01 1.0132e+00];  ps=ps/1013.25;
  end

pf=pt-ps;
ratio=pf./pt; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'/home/sergio/SPECTRA/PR_sigsig');

band   = 2350;
for ii = 1 : length(ps)
  for rr = 1 : 3
    [dsJJ(ii,rr),dfJJ(ii,rr),bsJJ(ii,rr),bfJJ(ii,rr),b2s,b2f] = ...
         co2_param_JJOHNS2002(band,pt(ii),ps(ii),rr);

    betaJJ(ii,rr) = (ps(ii)*bsJJ(ii,rr) + pf(ii)*bfJJ(ii,rr))/pt(ii);
    dJJ(ii,rr)    = (ps(ii)*dsJJ(ii,rr) + pf(ii)*dfJJ(ii,rr))/pt(ii);

    end
  end

fprintf(1,'   \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : length(ps)
  for rr = 1 : 3
    [dsRAL(ii,rr),dfRAL(ii,rr),bsRAL(ii,rr),bfRAL(ii,rr),b2s,b2f] = ...
       co2_param(band,pt(ii),ps(ii));

    betaRAL(ii,rr) = (ps(ii)*bsRAL(ii,rr) + pf(ii)*bfRAL(ii,rr))/pt(ii);
    dRAL(ii,rr)    = (ps(ii)*dsRAL(ii,rr) + pf(ii)*dfRAL(ii,rr))/pt(ii);

    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 1 : 3;    %there are 3 wavenumber regions

%loop over the pressure ratios
for ii = 1 : length(ps)
  plot(index,betaRAL(ii,:),'b+-',index,betaJJ(ii,:),'r+-');
  legend('RAL','JJ',0)
  tt = ['beta : ratio = ' num2str(ratio(ii))]; title(tt); pause; 
  end

index = 1 : 3;
for ii = 1 : length(ps)
  plot(index,dRAL(ii,:),'b+-',index,dJJ(ii,:),'r+-');
  legend('RAL','JJ',0)
  tt = ['doc : ratio = ' num2str(ratio(ii))]; title(tt); pause; 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wave1 = 2380:2391;  l1 = length(wave1);
wave2 = 2391:2400;  l2 = length(wave2);
wave3 = 2400:2405;  l3 = length(wave3);

wave = [wave1 wave2 wave3];
ii1 = (1 : l1);
ii2 = (1 : l2) + l1;
ii3 = (1 : l3) + l1 + l2;

for ii = 1 : length(ps)
  poRAL(ii1) = betaRAL(ii,1);
  poRAL(ii2) = betaRAL(ii,2);
  poRAL(ii3) = betaRAL(ii,3);
  poJJ(ii1) = betaJJ(ii,1);
  poJJ(ii2) = betaJJ(ii,2);
  poJJ(ii3) = betaJJ(ii,3);
  plot(wave,poRAL,'b+-',wave,poJJ,'r+-');
  legend('RAL','JJ',0)
  tt = ['beta : ratio = ' num2str(ratio(ii))]; title(tt); pause; 
  end

for ii = 1 : length(ps)
  poRAL(ii1) = dRAL(ii,1);
  poRAL(ii2) = dRAL(ii,2);
  poRAL(ii3) = dRAL(ii,3);
  poJJ(ii1) = dJJ(ii,1);
  poJJ(ii2) = dJJ(ii,2);
  poJJ(ii3) = dJJ(ii,3);
  plot(wave,poRAL,'b+-',wave,poJJ,'r+-');
  legend('RAL','JJ',0)
  tt = ['doc : ratio = ' num2str(ratio(ii))]; title(tt); pause; 
  end

