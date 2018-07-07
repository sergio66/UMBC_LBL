function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,layeramt,...
     prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO)

%[f,y]=yrun(temperature,band,ptotal,pself,layeramt,frequser,prb,LVF,IO,birn)
%this function computes CO2 R-line mixing due to the 2350 band
%disp('bands: 2350')

if ((band ~= 2350)&(band ~= 2351)&(band ~= 2352)&(band ~=2353)&(band~=2354))
  error('need 2350,2351,2352,2353,2354 band!!!!!!!!!')
  end

if (prb == 'p')
  prb='P';
  end
if (prb == 'r')
  prb='R';
  end

MGC=8.314674269981136; %%%%%%%%%%% correct value 
%layer amt in kmoles/cm^2
%density = n/V = partpress/(RT)
%density * path length = layer amt
%path length = layer amt/ density
path_length=layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2
den=101325*pself/MGC/temperature;    %density in no of moles/m^3
path_length=path_length/den*100;      %path length changed from m to cm

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie 
%dont forget to change pressure from torr to atm (p/760) 
%path_length=layeramt; 

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff]=...
    loader(temperature,band,path_length,ptotal,pself,prb);

fudger_and_threeregions;

elower        = elowerr;
[elower,jall] = efitter(jr,band,elowerr,elower,prb); 

[f,voi,fc,theRATIO] = ...  
    driver4um(temperature,ptotal,pself,layeramt,...   %gas layer params  
                  band,prb,LVF,IO,birn,...               %spectra computation  
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers  
                  path_length,elower,jall,...            %from initial call  
                  jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff); 
