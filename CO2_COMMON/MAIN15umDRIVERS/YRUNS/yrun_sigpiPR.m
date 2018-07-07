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
if ((band ~= 618)&(band ~= 648)&(band ~= 662)&...
    (band ~= 667)&(band ~= 720)&(band ~= 791)&(band ~= 2080))
  error('need bands to be 618,648,662,667,720,791,2080')
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
WOW=720;
[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb);

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

if (band == WOW)

  elower=elowerq; 
  [elower,jall]=efitter(jq,band,elowerq,elower,prb);  

  [W_co2for,W_co2self]=wfunco2er(jq,elower,elowerq,w_selfq,w_forq,band,...
                               jall,temperature,stuff); 
  [trans_ampl,population_t]=... 
        trans_pop(temperature,freqq,jq,elowerq,strenq,stuff); 

  W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref; 

  % multiply off diagonals of W by beta
  % this is new, and conforms to the code in
  %/salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run
  %and is done because of the last line in wfun1co2.m is
  %   dif=K+diag(wq_tape);     and not
  %   dif=beta*K1+diag(wq_tape);
  % this is to test no off diagonal elements vs lorentz
  % W_plus=diag(diag(W_plus));
  beta=stuff.beta;
  W_plus=W_plus.*(beta+(1-beta)*eye(length(jq))); 

  ratio=doratio(population_t,trans_ampl,W_plus); 
  if (ratio < 0)
    ratio=1e-10;         %if the sum rule gives -ve numbers
    end
  theRATIO=ratio;

  %now find out where the boundaries to do full mixing and just doing   
  %ratio to Lorentz, ie  
  [hinum,histart,histop,outnum,outstart,outstop,theratio] = ... 
    tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqq,ratio);  

  end

%%% only WOW band has linemixing computed in full detail
if (band == WOW)
  for kk=1:3
    theratio(kk)=ratio;

    if hinum(kk) > 0
    
      if (kk == 3)
        %if kk==3, do computations at high res, then boxcar integrate
        indHI = histart(kk):histop(kk);
        ind = outstart(kk):outstop(kk);
        ffff=f_hi(indHI);
      else
        %else computations at output res
        indHI = outstart(kk):outstop(kk);          %will NOT be used
        ind   = outstart(kk):outstop(kk);
        ffff=f_low(ind);
        end
      
      if ((LVF == 'F') | (LVF == 'f')) 
        if (kk == 3)          %do full mixing at high resolution
          voivoi=full2(freqq,ffff,W_plus,jq,w_selfq,w_forq,temperature,...
                trans_ampl,population_t,stuff,birn);
          scum=boxint2(voivoi,boxcar);  
          voi(ind)=scum;
        else                  %do mixing = lorentz*ratio at output resolution
          ymix=zeros(size(freqq)); 
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                 w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));   
          voi(ind)=voivoi;
          end

      elseif ((LVF == 'B') | (LVF == 'b'))  
        %do the fullmix first 
        if (kk == 3)         %do full mixing at high resolution
          voivoi=full2(freqq,ffff,W_plus,jq,w_selfq,w_forq,temperature,...
                  trans_ampl,population_t,stuff,birn);   
          voif=boxint2(voivoi,boxcar);  
        else                  %do mixing = lorentz*ratio at output res
          ymix=zeros(size(freqq)); 
          voif=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));   
          end     
        %do the first mix second 
        if (kk == 3)        %do first mixing at high reolution
          ymix=y1s(jq,freqq,elowerq,strenq,W_plus,trans_ampl,stuff,beta); 
          voivoi=voigtmix2(freqq,ffff,ymix,jq,...  
               temperature,w_forq,w_selfq,strenqt,stuff,layeramt,'V',birn);   
          voi1=boxint2(voivoi,boxcar);  
        else                  %do mixing = lorentz*ratio at output res
          ymix=zeros(size(freqq)); 
          voi1=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
              w_selfq,strenqt,stuff,layeramt,'V',birn,theratio(kk));   
          end
        voi(ind)=(ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

      else  
        %%%%first order line mixing calculation using voigt/lor lineshape  
        if ((IO == '1') & (kk == 3))       %do mixing at high resolution
          ymix=y1s(jq,freqq,elowerq,strenq,W_plus,trans_ampl,stuff,beta); 
          voivoi=voigtmix2(freqq,ffff,ymix,jq,...  
                temperature,w_forq,w_selfq,strenqt,stuff,layeramt,LVF,birn);   
          scum=boxint2(voivoi,boxcar);  
          voi(ind)=scum; 
        elseif ((IO == '1') & (kk ~= 3)) %mixing = lorentz*ratio at output res
          %zeroth order line mixing calculation using voigt/lor lineshape  
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                   w_selfq,strenqt,stuff,layeramt,LVF,birn,theratio(kk));   
          voi(ind)=voivoi;
        elseif ((IO == '0') & (kk == 3))     %do lorentz at high  res
          %zeroth order line mixing calculation using voigt/lor lineshape  
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                   w_selfq,strenqt,stuff,layeramt,LVF,birn,1.0);   
          scum=boxint2(voivoi,boxcar);  
          voi(ind)=scum;
        elseif ((IO == '0') & (kk ~= 3))     %do lorentz at output res
          %zeroth order line mixing calculation using voigt/lor lineshape  
          ymix=zeros(size(freqq));  
          voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
                w_selfq,strenqt,stuff,layeramt,LVF,birn,1.0);   
          voi(ind)=voivoi;
          end
        end         %if LVF = f or l or b
      end           %if hinum(ind) > 0
    end             %for kk=1:3
  %now smooth the transition from full mixing to RATIO mixing 
  %smooth over about 0.1 = 2 * 0.05 cm-1 
  df=abs(outwave(10)-outwave(1))/10; df=floor(0.05/df); 
  [voi]=smooth_full_ratio(voi,df,outnum,outstart,outstop); 

else
%%% other bands have "simple" linemixing
  %now find out where the boundaries to do full mixing and just doing   
  %ratio to Lorentz, ie  
  ratio=theRATIO;
  fprintf(1,'directly computing lineshape for PR band = %3i\n',band);
  fprintf(1,'if cousin not asked for, uses k/klor = %8.6f \n',theRATIO);
  [hinum,histart,histop,outnum,outstart,outstop,theratio] = ... 
    tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqq,ratio); 

%this is the code from doVmixSimple.f
%        if (nio .gt. 0) then !simple linemixing 
%          do j=1,lenf 
%            if (abs(v(j)-v0(i)) .gt. 15.0) then 
%              temp(j)=0.5*temp(j) 
%              endif 
%            end do 
%          endif 
 
  ymix=zeros(size(freqq)); 
  LVFtemp=LVF; 
  IOtemp=IO;
  if ((LVF == 'F') | (LVF == 'f')) 
    LVFtemp = 'V';
    IOtemp='1';
    end 
  if ((LVF == 'B') | (LVF == 'b')) 
    LVFtemp = 'V';
    IOtemp='1';
    end 
  for kk=1:3
    if hinum(kk) > 0
      if (kk == 3)
        %if kk==3, do computations at high res, then boxcar integrate
        indHI = histart(kk):histop(kk);
        ind = outstart(kk):outstop(kk);
        ffff=f_hi(indHI);
      else
        %else computations at output res
        indHI = outstart(kk):outstop(kk);          %will NOT be used
        ind   = outstart(kk):outstop(kk);
        ffff=f_low(ind);
        end
      voivoi=voigt_mix_simple2(freqq,ffff,ymix,jq,temperature,...
           w_forq,w_selfq,strenqt,stuff,layeramt,LVFtemp,IOtemp,birn,theRATIO);
      if (kk == 3)
        scum=boxint2(voivoi,boxcar);
        voi(ind)=scum;
      else
        voi(ind)=voivoi;
        end

      end             %if hinum(kk) > 0
    end               %for kk=1,3
  end                 %if band == WOW

fc=freqq;
cc=cputime-cc;
