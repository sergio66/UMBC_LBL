function [f,voi,fc,theRATIO]=yrun_sigpi662(temperature,band,...
  ptotal,pself,layeramt,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO) 
%function [f,mix4]=yrun662(temp,band,ptotal,pself,layeramt,frequser,NFM) 
%this function computes CO2 Q-line mixing due to the SigPi 662 bands  
%NFM is a latter variable that is N for NO line mixing
%                                 I for first order vanhub line mixing
%                                 F for full line mixing

fc = [];
cc=cputime;

MGC=8.314612138206606; %%%%%%%%%%% for testing against GLAB
MGC=8.314674269981136; %%%%%%%%%%% correct value
path_length=layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2 
den=101325*pself/MGC/temperature;%density in no of moles/m^3 
path_length=path_length/den*100;      %path length changed from m to cm 

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie
%dont forget to change pressure from torr to atm (p/760)
%path_length=layeramt; 

voi =zeros(size(outwave));
voi1=zeros(size(outwave));
voi2=zeros(size(outwave));

%if ((LVF == 'B') | (LVF == 'b')) 
%  voif=voi;
%  voi1=voi;
%  end

for ii=1:2
  [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader662(temperature,band,path_length,ptotal,pself,ii); 

  %oh boy! HAVE to do mixing computations
  elower=elowerq; 
  [elower,jall]=efitter(jq,band,elowerq,elower);  

  if ii==1
    [f1,voi1,fc1,theRATIO]=... 
      driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params  
                  band,'Q',LVF,IO,birn,...             %spectra computation  
                  boxcar,hiResFreq,outwave,theRATIO,...%boxcar, wavenumbers  
                  path_length,elower,jall,...          %from initial call  
                  jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,...
                  0,0,-1);                                %this is sigpi
%    [f1,voi1,fc1,theRATIO]=... 
%      driver15umPQR(temperature,ptotal,pself,layeramt,...  %gas layer params 
%                  band,'Q',LVF,IO,birn,...               %spectra computation 
%                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
%                  path_length,elower,jall,...            %from initial call 
%                  jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,...
%                  0,0,-1);                                %this is sigpi
    fc1=union(fc1,freqq);        

  elseif ii==2
    [f2,voi2,fc2,theRATIO]=... 
      driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params  
                  band,'Q',LVF,IO,birn,...             %spectra computation  
                  boxcar,hiResFreq,outwave,theRATIO,...%boxcar, wavenumbers  
                  path_length,elower,jall,...          %from initial call  
                  jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,...
                  0,0,-1);                                %this is sigpi 
%    [f2,voi2,fc2,theRATIO]=... 
%      driver15umPQR(temperature,ptotal,pself,layeramt,...  %gas layer params 
%                  band,'Q',LVF,IO,birn,...               %spectra computation 
%                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers 
%                  path_length,elower,jall,...            %from initial call 
%                  jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,...
%                  0,0,-1);                                %this is sigpi
    fc2=union(fc2,freqq);        

    end
  end

fc=union(fc1,fc2);
f=union(f1,f2);
voi=voi1+voi2;

