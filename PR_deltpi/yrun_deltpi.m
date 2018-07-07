function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,layeramt,...
                           prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO) 
%[f,mix4]=yrun(temperature,band,ptotal,pself,layeramt,freq,prb,LVF,IO,birn) 
%this function computes CO2 PR-line mixing due to the 668 and 740 bands  

%all this is put into structure "stuff" 
%global bsm pressure_self pressure_for btz B0 temperature_ref;   
%global density Boltzmann mass_CO2 speed_light path_length pressure_ref;  
%global K_scale_mixing K_scale_voigt K_scale_lor;  
%global population population_t t_rawdata voi_back voi_pr;  
%global beta_delt_self beta_delt_air beta_pi_self beta_pi_air;   

if ((band ~= 668)&(band ~= 740))
  error('need bands to be 668,740')
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

cc=cputime;

%%%WOW=-1111111111111111111;   %%%for PR_deltpi, band to do full mixing!!!
%%%for PR_deltpi, originally had NO band on which we did line mixing; just
%%%used k/klor=0.5 or whatever fav number.

%%%but we notice things get bad on the 740 part of the spectrum
WOW=740;     %%%%%%for PR_deltpi, band to do full mixing!!!

[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb);

if (band == WOW)
  elower=elowerq; 
  [elower,jall]=efitter(jq,band,elowerq,elower,prb);  
  beta_pi=(pself*stuff.beta_pi_self+(ptotal-pself)*stuff.beta_pi_air)/ptotal;
  beta_delt=(pself*stuff.beta_delt_self+...
                (ptotal-pself)*stuff.beta_delt_air)/ptotal;
  [f,voi,fc,theRATIO]=... 
     driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params 
                band,prb,LVF,IO,birn,...               %spectra computation 
                boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
                path_length,elower,jall,...            %from initial call 
                jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,... 
                beta_pi,beta_delt,1);    %sigpi,deltpi slightly differ 

else
%%%     driver15umPR_SIMPLE(temperature,ptotal,pself,layeramt,... %gas params 
  [f,voi,fc,theRATIO]=... 
    driver15umPR(temperature,ptotal,pself,layeramt,... %gas params 
                band,prb,LVF,IO,birn,...               %spectra computation 
                boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
                path_length,...            %from initial call 
                jq,w_forq,w_selfq,freqq,strenqt,stuff);
  end

