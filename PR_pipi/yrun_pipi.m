function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,...
               layeramt,prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO)
%func [f,fullmix4]=yrun(temperature,band,ptotal,pself,layeramt,frequser,prb)
%this function computes CO2 R-line mixing due to the 2320 band
%disp('bands: 2320')

if ((band ~= 2320) & (band ~= 2321) & (band ~= 2322))
  error('need 2320,2321, 2322 band!!!!!!!!!')
  end

if (prb == 'p')
  prb = 'P';
  end
if (prb == 'r')
  prb = 'R';
  end

MGC = 8.314674269981136; %%%%%%%%%%% correct value 
% layer amt in kmoles/cm^2
% density  =  n/V  =  partpress/(RT)
% density * path length  =  layer amt
% path length  =  layer amt/ density
path_length = layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2
den = 101325*pself/MGC/temperature;     %density in no of moles/m^3
path_length = path_length/den*100;      %path length changed from m to cm

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie 
%dont forget to change pressure from torr to atm (p/760) 
%path_length = layeramt; 

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff] = ...
    loader(temperature,band,path_length,ptotal,pself,prb);

fudger_and_threeregions;

elower = elowerr;
if ((band  == 2320) | (band == 2321))
  [elower,jall] = efitter2320(jr,band,elowerr,elower,prb); 
elseif (band  == 2322)
  [elower,jall] = efitter2322(jr,band,elowerr,elower,prb); 
  end

[f,voi,fc,theRATIO] = ...   
  driver4um(temperature,ptotal,pself,layeramt,...   %gas layer params   
            band,prb,LVF,IO,birn,...               %spectra computation   
            boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers   
            path_length,elower,jall,...            %from initial call   
           jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff);  
 
