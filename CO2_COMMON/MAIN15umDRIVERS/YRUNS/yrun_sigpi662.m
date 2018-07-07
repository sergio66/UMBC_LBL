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

voi=zeros(size(outwave));
if ((LVF == 'B') | (LVF == 'b')) 
  voif=voi;
  voi1=voi;
  end

for ii=1:2
  [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader662(temperature,band,path_length,ptotal,pself,ii); 

  f_hi=hiResFreq;
  f_low=outwave;
  f=outwave;

  elower=elowerq; 

  %oh boy! HAVE to do mixing computations
  [elower,jall]=efitter(jq,band,elowerq,elower);  
  [W_co2for,W_co2self]=wfunco2er(jq,elower,elowerq,w_selfq,w_forq,band,...
                                 jall,temperature,stuff); 
  [trans_ampl,population_t]=... 
        trans_pop(temperature,freqq,jq,elowerq,strenq,stuff); 

  W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref; 
  
  beta=stuff.beta;
  W_plus=W_plus.*(beta+(1-beta)*eye(length(jq))); 

  ratio=doratio(population_t,trans_ampl,W_plus); 
  if (ratio <= 0.0)  
    ratio=1.0e-10;         %if the sum rule give -ve numbers 
    end
  theRATIO=ratio;

  %now find out where the boundaries to do full mixing and just doing  
  %ratio to Lorentz, ie 
  [hinum,histart,histop,outnum,outstart,outstop,theratio] = ...
    tempratio15um(band,temperature,boxcar,f_hi,f_low,'Q',freqq,ratio);

  for kk=1:3
    theratio(kk)=ratio;

    if (hinum(kk) > 0)

      if (kk == 3)
        %if kk==3, do computations at high res, then boxcar integrate
        indHI = histart(kk):histop(kk);
        ind   = outstart(kk):outstop(kk);
        ffff=f_hi(indHI);
      else
        %else computations at output res
        indHI = outstart(kk):outstop(kk);          %will NOT be used
        ind   = outstart(kk):outstop(kk);
        ffff=f_low(ind);
        end

      if ((LVF == 'F') | (LVF == 'f')) 
        if (kk==3)
          voivoi=full2(freqq,ffff,W_plus,jq,w_selfq,w_forq,temperature,...
              trans_ampl,population_t,stuff,birn);  
          scum=boxint2(voivoi,boxcar);
          voi(ind)=voi(ind)+scum;
        else                  %do mixing = lorentz*ratio 
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
               w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));    
          voi(ind)=voi(ind)+voivoi; 
          end 
 
      elseif ((LVF == 'B') | (LVF == 'b'))  
        %do the fullmix first 
        if (kk==3)
          voivoi=full2(freqq,ffff,W_plus,jq,w_selfq,w_forq,temperature,...
                trans_ampl,population_t,stuff,birn);   
          scum=boxint2(voivoi,boxcar);  
          voif(ind)=voif(ind)+scum;
        else                  %do mixing = lorentz*ratio 
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
               w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));    
          voif(ind)=voif(ind)+voivoi;
          end      
        %do the first mix second 
        if (kk==3)           %do fist mixing
          ymix=y1s(jq,freqq,elowerq,strenq,W_plus,trans_ampl,stuff,beta); 
          voivoi=voigtmix2(freqq,ffff,ymix,jq,...  
               temperature,w_forq,w_selfq,strenqt,stuff,layeramt,'V',birn);   
          scum=boxint2(voivoi,boxcar);  
          voi1(ind)=voi1(ind)+scum;
        else 
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
              w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));    
          voi1(ind)=voi1(ind)+voivoi; 
          end 
         voi(ind)=voi(ind) + voi1(ind) + ... 
            (ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif(ind)-voi1(ind));

      else  
        %%%%first order line mixing calculation using voigt/lor lineshape 
        if ((IO == '1') & (kk == 3))
          ymix=y1s(jq,freqq,elowerq,strenq,W_plus,trans_ampl,stuff,beta);
          voivoi=voigtmix2(freqq,ffff,ymix,jq,... 
                 temperature,w_forq,w_selfq,strenqt,stuff,layeramt,LVF,birn);  
          scum=boxint2(voivoi,boxcar);  
          voi(ind)=voi(ind)+scum;
        elseif ((IO == '1') & (kk ~= 3))
          %zeroth order line mixing calculation using voigt/lor lineshape 
          ymix=zeros(size(freqq)); 
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
              w_selfq,strenqt,stuff,layeramt,LVF,birn,theratio(kk));  
          voi(ind)=voi(ind)+voivoi;
        elseif ((IO == '0') & (kk == 3))
          %zeroth order line mixing calculation using voigt/lor lineshape 
          ymix=zeros(size(freqq)); 
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                   w_selfq,strenqt,stuff,layeramt,LVF,birn,1.00);  
          scum=boxint2(voivoi,boxcar);  
          voi(ind)=voi(ind)+scum;
        elseif ((IO == '0') & (kk ~= 3))
          %zeroth order line mixing calculation using voigt/lor lineshape 
          ymix=zeros(size(freqq)); 
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
              w_selfq,strenqt,stuff,layeramt,LVF,birn,1.00);  
          voi(ind)=voi(ind)+voivoi;
          end
 
       end           %if LVF=b,v,f
     fc=union(fc,freqq);
     end             %if hinum(kk) > 0
   end               %for kk=1,3
  %now smooth the transition from full mixing to RATIO mixing 
  %smooth over about 0.1 cm-1 
  df=abs(outwave(10)-outwave(1))/10; df=floor(0.05/df); 
  [voi]=smooth_full_ratio(voi,df,outnum,outstart,outstop); 
  end                 %for ii=1,2

cc=cputime-cc;
