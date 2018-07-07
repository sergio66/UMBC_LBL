function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,...
            layeramt,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO)
%function [f,fullmix4]=yrun(temperature,band,ptotal,pself,layeramt,frequsr,NIF)
%this function computes CO2 Q-line mixing due to the 668,740,2093 bands 
%disp('bands: 668, 740, 2093')

% boxcar,hiResFreq,outwave  ... are to do all the boxcar integrations HERE!!!!

%all this is put into structure "stuff"
%global bsm pressure_self pressure_for btz B0 temperature_ref;  
%global density Boltzmann mass_CO2 speed_light path_length pressure_ref; 
%global K_scale_mixing K_scale_voigt K_scale_lor; 
%global population population_t t_rawdata voi_back voi_pr; 
%global beta_delt_self beta_delt_air beta_pi_self beta_pi_air;  

if ((band ~= 668) & (band ~= 740) & (band~=2093))
  error('need bands to be 668, 740 or 2093');
  end

%http://physics.nist.gov/cuu/Constants/index.html?/codata86.html
%note from NIST : 8.314510=1.38e-23*6.023e23
%thus avogadro = 6.0221367e23, boltzmann = 1.380658e-23
%MGC= 8.314510;

%but these are rather inconsistent with Dave Edwards gendat.f
%  /salsify/packages/Genln2/Genln2/gendat.f
%       DATA TS/296.0/ PS/1.0/
%       DATA 
%     + C1/1.1911E-8/ C2/1.4387863/ VLIGHT/2.9979245E8/ 
%     + AVOG/6.022045E26/ R2/11526.218/ GRAV/9.80665/
%     + CPAIR/1006/ ATMB/1013.25/
%       DATA PI/3.1415926/ 
%we know c2=hc/kb 
%we know c1=2hc*c ==> h =  6.626387762479327e-26
%thus kb that dave edwards uses = 1.380706100665328e-23
%thus MGC should be 8.314674269981136
%%%%%%%%%%%% but in reading in HITLIN file, he multiplies line strength by
%%%%%%%%%%%% 6.022e26 instead of avog!!!!!!!!!!
% ie he uses 9.999925274553744e-01 times this amount

cc=cputime;

MGC=8.314612138206606; %%%%%%%%%%% for testing against GLAB
MGC=8.314674269981136; %%%%%%%%%%% correct value

%layer amt in kmoles/cm^2
%density = n/V = partpress/(RT)
%density * path length = layer amt
%path length = layer amt/ density
path_length = layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2
den         = 101325*pself/MGC/temperature; %density in no of moles/m^3
path_length = path_length/den*100;          %path length changed from m to cm

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie
%dont forget to change pressure from torr to atm (p/760)
%path_length=layeramt;

[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=...
    loader(temperature,band,path_length,ptotal,pself);

elower        = elowerq;
[elower,jall] = efitter(jq,band,elowerq,elower); 

beta_pi   = (pself*stuff.beta_pi_self   + (ptotal-pself)*stuff.beta_pi_air);
beta_pi   = beta_pi/ptotal;
beta_delt = (pself*stuff.beta_delt_self + (ptotal-pself)*stuff.beta_delt_air);
beta_delt = beta_delt/ptotal;

[f,voi,fc,theRATIO]=... 
      driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params 
                 band,'Q',LVF,IO,birn,...               %spectra computation 
                 boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
                 path_length,elower,jall,...            %from initial call 
                 jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,... 
                 beta_pi,beta_delt,1);   %sigpi,deltpi slightly differ 
fc=freqq;
cc=cputime-cc;
