function [f,voi,fc,theRATIO]=...
        driver15umPQR(temperature,ptotal,pself,layeramt,...  %gas layer params
                  band,prb,LVF,IO,birn,...               %spectra computation
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers
                  path_length,elower,jall,...            %from initial call
                  jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff,...
                  beta_pi,beta_delt,sigORdelt)    %sigpi,deltpi slightly differ

%this is for the 700cm-1  (15 um band)
%
%  ----------------|-----------------------------------------|------------- 
%                  f1                                        f2 
%   k=klor*ratio                 k=kfull                        k=klor*ratio 
% 
%         I                        III                            II 
%
  
if (abs(sigORdelt) ~=1) 
  sigORdelt
  error('this driver is for the PQR sigpi,deltpi mixing calcs!!')
  end
  

%this section of code is identical for PQQ bands of the 15 um mixing code
cc=cputime;

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

if (sigORdelt < 0)              %this is sigpi
  [W_co2for,W_co2self]=wfunco2er(jq,elower,elowerq,w_selfq,w_forq,...
                                   band,jall,temperature,stuff);
  [trans_ampl,population_t]=...
      trans_pop(temperature,freqq,jq,elowerq,strenq,stuff);

  W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

  % multiply off diagonals of W by beta 
  % this is new, and conforms to the code in 
  % /salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run 
  %and is done because of the last line in wfun1co2.m is 
  %   dif=K+diag(wq_tape);     and not 
  %   dif=beta*K1+diag(wq_tape);   
  beta=stuff.beta; 
  W_plus=W_plus.*(beta+(1-beta)*eye(length(jq)));   
  %%% to test linewidths  W_plus=diag(diag(W_plus)); 

elseif (sigORdelt > 0)              %this is deltpi
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
  end

ratio=doratio(population_t,trans_ampl,W_plus);
if (ratio <= 0.0) 
  ratio=1.0e-10;         %if the sum rule give -ve numbers
  end
theRATIO=ratio;

%now find out where the boundaries to do full mixing and just doing    
%ratio to Lorentz, ie   
[hinum,histart,histop,outnum,outstart,outstop,theratio] = ...  
  tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqq,ratio);  
  
voi=zeros(size(outwave));

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
