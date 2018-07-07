ido = input('do the RAL (+1) or JJOHNS (-1) p/ps parameters : ');

%%%region == 1 ==> 2380-2391
%%%region == 2 ==> 2391-2400
%%%region == 3 ==> 2400-2405
region = input('Enter region (1,2,3) : ');

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
  [dsJJ(ii),dfJJ(ii),bsJJ(ii),bfJJ(ii),b2s,b2f] = ...
       co2_param_JJOHNS2002(band,pt(ii),ps(ii),region);
  end
betaJJ = (ps.*bsJJ + pf.*bfJJ)./pt;
dJJ    = (ps.*dsJJ + pf.*dfJJ)./pt;

fprintf(1,'   \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : length(ps)
  [dsRAL(ii),dfRAL(ii),bsRAL(ii),bfRAL(ii),b2s,b2f] = ...
       co2_param(band,pt(ii),ps(ii));
  end
betaRAL = (ps.*bsRAL + pf.*bfRAL)./pt;
dRAL    = (ps.*dsRAL + pf.*dfRAL)./pt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(ratio,dsRAL,'+-',ratio,dsJJ,'r+-')
tt = ['ds : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

plot(ratio,dfRAL,'+-',ratio,dfJJ,'r+-')
tt = ['df : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

plot(ratio,dRAL,'+-',ratio,dJJ,'r+-')
tt = ['doc : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

plot(ratio,bsRAL,'+-',ratio,bsJJ,'r+-')
tt = ['bs : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

plot(ratio,bfRAL,'+-',ratio,bfJJ,'r+-')
tt = ['bf : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

plot(ratio,betaRAL,'+-',ratio,betaJJ,'r+-')
tt = ['beta : Region ' num2str(region)];
title(tt); legend('ral','jjohns',0)
pause

[ratio; betaRAL./betaJJ; dRAL./dJJ]'