function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,layeramt,...
                         prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO) 
%[f,mix4]=yrun(temperature,band,ptotal,pself,layeramt,freq,prb,LVF,IO,birn) 
%this function computes CO2 PR-line mixing due to the 667 and 720 bands  

%all this is put into structure "stuff" 
%global bsm pressure_self pressure_for btz B0 temperature_ref;   
%global density Boltzmann mass_CO2 speed_light path_length pressure_ref;  
%global K_scale_mixing K_scale_voigt K_scale_lor;  
%global population population_t t_rawdata voi_back voi_pr;  
%global beta_delt_self beta_delt_air beta_pi_self beta_pi_air;   

%load /beet/users/tobin/Hittomat/Src/hit_co2_15um
PR_sigpibands=[618 648 662 667 720 791 1932 2080 2129];
len=length(intersect(band,PR_sigpibands));
if (len < 1)
  error('need bands to be 618,648,662,667,720,791,1932,2080,2129')
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

%%%%%%%%%%%%%%%%%%%WOW=-11111111111;
WOW1=720;    %%%%%%%this is the one where we do mixing; else use k/klor=const
WOW2=667;    %%%%%%%this is the one where we do mixing; else use k/klor=const
%%%% WOW2=720;    %%%%%%%this is the one where we do mixing; else use k/klor=const commented this out 6/16/2015
WOW3=791;    %%%%%%%this is the one where we do mixing; else use k/klor=const
[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb);

%%%%if ((band == WOW1) | (band==WOW2))
%%if ((band == WOW1)|(band==WOW2)|(band == WOW3))
if ((band == WOW1))
  elower=elowerq; 
  [elower,jall]=efitter(jq,band,elowerq,elower,prb);  
  [f,voi,fc,theRATIO]=... 
     driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params 
                   band,prb,LVF,IO,birn,...               %spectra computation 
                   boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
                   path_length,elower,jall,...            %from initial call 
                   jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,... 
                   0,0,-1);    %sigpi,deltpi slightly differ 

else
%%%%     driver15umPR_SIMPLE(temperature,ptotal,pself,layeramt,... %gas params
  [f,voi,fc,theRATIO]=... 
    driver15umPR(temperature,ptotal,pself,layeramt,...  %gas params 
                 band,prb,LVF,IO,birn,...               %spectra computation 
                 boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
                 path_length,...            %from initial call 
                 jq,w_forq,w_selfq,freqq,strenqt,stuff);
  end
