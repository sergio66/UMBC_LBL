function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,layeramt,...
     prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO)
%[f,y]=yrun(temperature,band,ptotal,pself,layeramt,frequser,prb,LVF,IO,birn)
%this function computes CO2 R-line mixing due to the 2350 band
%disp('bands: 2350')

%all this is put into structure "stuff"
%global bsm pressure_self pressure_for btz B0 temperature_ref;  
%global density Boltzmann mass_CO2 speed_light path_length pressure_ref; 
%global K_scale_mixing K_scale_voigt K_scale_lor; 
%global population population_t t_rawdata voi_back voi_pr; 
%global beta_delt_self beta_delt_air beta_pi_self beta_pi_air;  

if ((band ~= 2350)&(band ~= 2351)&(band ~= 2352)&(band ~=2353)&(band~=2354))
  error('need 2350,2351,2352,2353,2354 band!!!!!!!!!')
  end

if (prb == 'p')
  prb='P';
  end
if (prb == 'r')
  prb='R';
  end

cc=cputime;

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

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

elower=elowerr;

[elower,jall]=efitter(jr,band,elowerr,elower,prb); 

%beta=(pself*stuff.beta_pure+(ptotal-pself)*stuff.beta_for)./ptotal;
[W_co2for,W_co2self]=wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuff);
[trans_ampl,population_t]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuff);

W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

% multiply off diagonals of W by beta 
% this is new, and conforms to the code in 
%/salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run 
%and is done because of the last line in wfun1co2.m is 
%dif=K+diag(wq_tape);  and not 
%dif=beta*K1+diag(wq_tape); 
beta=stuff.beta;
W_plus=W_plus.*(beta+(1-beta)*eye(length(jr)));

%[lor,lormix,lor_birn,lormix_birn]=...
%klormix_birn(f,freqr,W_plus,jr,temperature,beta,w_selfr,w_forr,...
%trans_ampl,population_t,stuff,strenrt,ymix);
%fullmix4=full_mix4(f,freqr,W_plus,jr,temperature,beta,w_selfr,w_forr,...
%trans_ampl,population_t,stuff);

ratio=doratio(population_t,trans_ampl,W_plus);
if (ratio <= 0.0) 
  ratio=1.0e-10;         %if the sum rule give -ve numbers
  end
%%ZZ
%if (band ~=2350)
%  ratio=0.1;
%  end
theRATIO=ratio;

%now find out where the boundaries to do full mixing and just doing 
%ratio to Lorentz, ie
[hinum,histart,histop,outnum,outstart,outstop,theratio]=...
   tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqr,ratio);

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
        voivoi=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
             trans_ampl,population_t,stuff,birn);
        scum=boxint2(voivoi,boxcar);  
        voi(ind)=scum;
      else                  %do mixing = lorentz*ratio at output resolution
        ymix=zeros(size(freqr)); 
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,...
            w_selfr,strenrt,stuff,layeramt,'V',birn,theratio(kk));
        voi(ind)=voivoi;
        end

    elseif ((LVF == 'B') | (LVF == 'b'))  
      %do the fullmix first 
      if (kk == 3)         %do full mixing at high resolution
        voivoi=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
              trans_ampl,population_t,stuff,birn);   
        voif=boxint2(voivoi,boxcar);  
      else                  %do mixing = lorentz*ratio at output res
        ymix=zeros(size(freqr)); 
        voif=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn,theratio(kk));
        end 

      %do the first mix second 
      if (kk == 3)        %do first mixing at high reolution
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);   
        voi1=boxint2(voivoi,boxcar);  
      else                  %do mixing = lorentz*ratio at output res
        ymix=zeros(size(freqr)); 
        voi1=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                 w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn,theratio(kk));

        end
      voi(ind)=(ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

    else  
      %%%%first order line mixing calculation using voigt/lor lineshape  
      %%%   for first order mixing, usingvoigtmixRatioQUAD the last 3 
      %%%   parameters would be  theratio(kk),theslope(kk),theendpt(kk));
      %%%   for NO mixing, usingvoigtmixRatioQUAD the last 3 
      %%%   parameters would be  1.0,0.0,theendpt(kk));      
      if ((IO == '1') & (kk == 3))       %do mixing at high resolution
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);   
        scum=boxint2(voivoi,boxcar);  
        voi(ind)=scum; 
      elseif ((IO == '1') & (kk ~= 3))%do mixing = lorentz*ratio at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,theratio(kk));
        voi(ind)=voivoi;
      elseif ((IO == '0') & (kk == 3))     %do lorentz at high  res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                 w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,1.0);
        scum=boxint2(voivoi,boxcar);  
        voi(ind)=scum;
      elseif ((IO == '0') & (kk ~= 3))     %do lorentz at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,1.0);
        voi(ind)=voivoi;
        end
      end         %if LVF = f or l or b
    end           %if hinum(ind) > 0
  end


%now smooth the transition from full mixing to RATIO mixing
%smooth over about 0.1 cm-1
df=abs(outwave(10)-outwave(1))/10; df=floor(0.05/df);
[voi]=smooth_full_ratio(voi,df,outnum,outstart,outstop);

%check to make sure none of the terms are negative ....
ind=1:50:length(voi);
voiind=voi(ind);
b=sum(voiind<0); %if any of the terms are negative, have to fix
dave=-1;
if ((b > 0) & (dave < 0))
  fprintf(1,'some of the k''s are negative .. resetting to 0 \n');
  clear b
  b=find(voi < 0);
  voi(b) = 0.0;
  end

fc=freqr;
cc=cputime-cc;

