function [thediff] = fit_run7(X,topts,labdata,sim)

%%% X = lsqnonlin(@fit_run7,X,[],[],sopts,topts,labdata,sim);

global betafudge docfudge

fprintf(1,'in fit_run7, X = [%8.6f %8.6f %8.6f ]\n',X(1),X(2),X(3));

betafudge = X(1);
docfudge  = X(2);
freqshift = X(3);
f         = labdata.f;
t_rawdata = labdata.data;
backgnd   = sim.backgnd;
pipi_deltdelt_R = sim.pp_dd_KC;  %%%kcarta simulation at 0.0025 spacing 
kweakco2        = sim.weakco2;   %%%kcarta simulation at 0.0025 spacing 
kn2             = sim.n2_KC;     %%%kcarta simulation at 0.0025 spacing 

ipfileIN  = sim.ipfile;
if (sum(isnan(X)) > 0)
  thediff = ones(length(f));
elseif (X(1) < 0.02)
  thediff = ones(length(f));
elseif (X(1) > 2)
  thediff = ones(length(f));
elseif (X(2) < 0.02)
  thediff = ones(length(f));
elseif (X(2) > 5)
  thediff = ones(length(f));
elseif (abs(X(3)) > 1)
  thediff = ones(length(f));
else 

  path(path,'/home/sergio/SPECTRA');
  cd /home/sergio/SPECTRA
  [fr,kSS] = run7co2_JJOHNS(2,2255,2305,ipfileIN,topts);
  cd /home/sergio/SPECTRA/CO2_COUSIN_2351
  k = interp1(fr + X(3),kSS + pipi_deltdelt_R + kn2 + kweakco2,f);

  whos k backgnd t_rawdata
  totalT = exp(-(backgnd + k));
  thediff0 = t_rawdata - totalT;

  %%right now we wanna weight stuff equally
  iweight = +1;
  if iweight > 0
    weight = ones(size(thediff0));
    thediff = thediff0 .* weight;
    end

  %%right now we wanna weight stuff from 2280 - 2285 cm-1
  iweight = +1;
  if iweight > 0
    weight = ones(size(thediff0));
    popo = find((f >= 2275) & (f <= 2290));
    weight(popo) = 10.0;
    thediff = thediff0 .* weight;
    end

  plot(f,t_rawdata,f,totalT,f,thediff0,f,thediff);
  pause(0.1);  
  end