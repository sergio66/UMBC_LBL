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
path_length=layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2
den=101325*pself/MGC/temperature;    %density in no of moles/m^3
path_length=path_length/den*100;      %path length changed from m to cm

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie
%dont forget to change pressure from torr to atm (p/760)
%path_length=layeramt;

[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=...
    loader(temperature,band,path_length,ptotal,pself);

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

elower=elowerq;

%oh boy! HAVE to do mixing computations
[elower,jall]=efitter(jq,band,elowerq,elower); 
beta_pi=(pself*stuff.beta_pi_self+(ptotal-pself)*stuff.beta_pi_air)/ptotal;
beta_delt=(pself*stuff.beta_delt_self+...
                (ptotal-pself)*stuff.beta_delt_air)/ptotal;
[W_co2for,W_co2self]=wfunco2er(jq,elower,elowerq,w_selfq,w_forq,band,...
             jall,temperature,stuff,beta_pi,beta_delt);
[trans_ampl,population_t]=...
      trans_pop(temperature,freqq,jq,elowerq,strenq,stuff);

%disp('computing 1st order mixing coefficients for foriegn broadening')
%y_for=y1s(jq,freqq,elowerq,strenq,W_co2for,trans_ampl,beta_pi,beta_delt);
%disp('computing 1st order mixing coefficients for self broadening')
%y_self=y1s(jq,freqq,elowerq,strenq,W_co2self,trans_ampl,beta_pi,beta_delt);
%ymix=(pself*y_self+(ptotal-pself)*y_for)/stuff.pressure_ref;
%[lor,lormix]=klormix(freqq,f,ymix,jq,temperature,w_forq,w_selfq,...
%                     strenqt,stuff);

W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

%%multiply off diagonals by beta; this code has been removed from
%% fullmix4.m and inserted here

no_lines=length(freqq); 
inde=find(rem(jq,2)==0);indo=find(rem(jq,2)==1); 
x=ones(no_lines)-eye(no_lines); 
x(inde,inde)=x(inde,inde)*abs(beta_pi*beta_delt); 
x(indo,indo)=x(indo,indo)*abs(beta_pi*beta_delt); 
x(inde,indo)=x(inde,indo)*abs((1-beta_pi)*(1-beta_delt)); 
x(indo,inde)=x(indo,inde)*abs((1-beta_pi)*(1-beta_delt)); 
x=x+eye(no_lines); 
W_plus=W_plus.*x; 
beta=beta_pi;  %just a dummy for y1s

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

  if hinum(kk) > 0
    
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
        voif=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,w_selfq,...
                strenqt,stuff,layeramt,'V',birn,theratio(kk));   
        end     
      %do the first mix second 
      if (kk == 3)        %do first mixing at high reolution
        ymix=y1s(jq,freqq,elowerq,strenq,W_plus,trans_ampl,stuff,beta); 
        voivoi=voigtmix2(freqq,ffff,ymix,jq,...  
               temperature,w_forq,w_selfq,strenqt,stuff,layeramt,'V',birn);   
        voi1=boxint2(voivoi,boxcar);  
      else                  %do mixing = lorentz*ratio at output res
        ymix=zeros(size(freqq)); 
        voi1=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,w_selfq,...
                 strenqt,stuff,layeramt,'V',birn,theratio(kk));   
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
      elseif ((IO == '1') & (kk ~= 3))%do mixing = lorentz*ratio at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqq));  
        voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,w_selfq,...
              strenqt,stuff,layeramt,LVF,birn,theratio(kk));   
        voi(ind)=voivoi;
      elseif ((IO == '0') & (kk == 3))     %do lorentz at high  res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqq));  
        voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,w_selfq,...
                 strenqt,stuff,layeramt,LVF,birn,1.0);   
        scum=boxint2(voivoi,boxcar);  
        voi(ind)=scum;
      elseif ((IO == '0') & (kk ~= 3))     %do lorentz at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqq));  
        voivoi=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,w_selfq,...
              strenqt,stuff,layeramt,LVF,birn,1.0);   
        voi(ind)=voivoi;
        end
      end         %if LVF = f or l or b
    end           %if hinum(ind) > 0
  end             %for kk=1:3

%now smooth the transition from full mixing to RATIO mixing 
%smooth over about 0.1 cm-1 
df=abs(outwave(10)-outwave(1))/10; df=floor(0.05/df); 
[voi]=smooth_full_ratio(voi,df,outnum,outstart,outstop); 
  
fc=freqq;
cc=cputime-cc;
