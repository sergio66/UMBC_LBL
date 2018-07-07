function [f,voi,fc,theRATIO]=...
        driver15um(temperature,ptotal,pself,layeramt,...  %gas layer params
                  band,prb,LVF,IO,birn,...               %spectra computation
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers
                  path_length,elower,jall,...            %from initial call
                  jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff,...
                  beta_pi,beta_delt,sigORdelt)    %sigpi,deltpi slightly differ

%this is really the same as driver4um.m except that it allows in 3 additional 
%parameters, so that it can do the Q_deltpi linemixing plus birnbaum

%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%

global quiet 

if (abs(sigORdelt) ~=1) 
  sigORdelt
  error('this driver is for the PQR sigpi,deltpi mixing calcs!!')
  end
  
%this section of code is identical for PQQ bands of the 15 um mixing code
cc=cputime;

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

stuff.temperature = temperature;
stuff.prb         = prb;

if (sigORdelt < 0)              %this is sigpi
  [W_co2for,W_co2self]=wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuff);
  [trans_ampl,population_t]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuff);

  W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

  % multiply off diagonals of W by beta 
  % this is new, and conforms to the code in 
  % /salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run 
  %and is done because of the last line in wfun1co2.m is 
  %   dif=K+diag(wq_tape);     and not 
  %   dif=beta*K1+diag(wq_tape);   
  beta=stuff.beta; 
  W_plus=W_plus.*(beta+(1-beta)*eye(length(jr)));   
  %%% to test linewidths  W_plus=diag(diag(W_plus)); 

elseif (sigORdelt > 0)              %this is deltpi
  [W_co2for,W_co2self]=wfunco2er(jr,elower,elowerr,w_selfr,w_forr,band,... 
             jall,temperature,stuff,beta_pi,beta_delt); 
  [trans_ampl,population_t]=... 
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuff); 
 
  %disp('computing 1st order mixing coefficients for foriegn broadening') 
  %y_for=y1s(jr,freqr,elowerr,strenr,W_co2for,trans_ampl,beta_pi,beta_delt); 
  %disp('computing 1st order mixing coefficients for self broadening') 
  %y_self=y1s(jr,freqr,elowerr,strenr,W_co2self,trans_ampl,beta_pi,beta_delt); 
  %ymix=(pself*y_self+(ptotal-pself)*y_for)/stuff.pressure_ref; 
  %[lor,lormix]=klormix(freqr,f,ymix,jr,temperature,w_forr,w_selfr,... 
  %                     strenrt,stuff); 
 
  W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref; 
 
  %%multiply off diagonals by beta; this code has been removed from 
  %% fullmix4.m and inserted here 
 
  no_lines=length(freqr);  
  inde=find(rem(jr,2)==0);indo=find(rem(jr,2)==1);  
  x=ones(no_lines)-eye(no_lines);  
  x(inde,inde)=x(inde,inde)*abs(beta_pi*beta_delt);  
  x(indo,indo)=x(indo,indo)*abs(beta_pi*beta_delt);  
  x(inde,indo)=x(inde,indo)*abs((1-beta_pi)*(1-beta_delt));  
  x(indo,inde)=x(indo,inde)*abs((1-beta_pi)*(1-beta_delt));  
  x=x+eye(no_lines);  
  W_plus=W_plus.*x;  
  beta=beta_pi;  %just a dummy for y1s 
  end

if quiet > 0 
  fprintf(1,'beta= %8.6f d of c = %9.5f \n',beta,stuff.duration); 
  end 

ratio=doratio(population_t,trans_ampl,W_plus);
if (ratio <= 0.0) 
  ratio=1.0e-10;         %if the sum rule give -ve numbers
  end
theRATIO=ratio;

%now find out where the boundaries to do full mixing and just doing    
%ratio to Lorentz, ie   
[hinum,histart,histop,outnum,outstart,outstop,theratio] = ...  
  tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqr,ratio,...
              strenrt,w_forr,w_selfr,jr,stuff);

voi=zeros(size(outwave));

for kk=1:5
  theratio(kk)=ratio;

  if hinum(kk) > 0

    ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);

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
      if ((kk==1)|(kk==2)) 
        %do cousin if necessary
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          tempbirn = 'c';
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          end
        %note we use k/kfar=1.0 instead of theratio(kk)
        ymix=zeros(size(freqr)); 
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,...
            w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        voi(ind)=voivoi;

      elseif (kk == 3)          %do full mixing at high resolution
        scum=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                   trans_ampl,population_t,stuff,birn,theratio(kk),...
                   strenrt,layeramt,ymix,band,prb);
        voivoi=boxint2(scum,boxcar);  
        voi(ind)=voivoi;

      elseif ((kk==4)|(kk==5))      %do cousin + linemixing blend
        if quiet > 0
          fprintf(1,'doing cousin + linemix blend .. ratio = 1.0 \n');
          end
        voivoi=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb); 
        ymix=zeros(size(freqr)); 
        %note we use k/kfar=1.0 instead of theratio(kk)
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          tempbirn = 'c';
          end
        voivoiC=voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,...
            w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        %these are the start and stop points of current wavevector
        s1=ffff(1);   s2=ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend]=blender15um(band,prb,freqr,kk);
        if (kk==4)       %%%%%left side
          %%%%%%%%%%%%%%   COUS       |     BLEND        |       FULL
          %%%%%%%%%%%%%               fa                fb
          fb=endfull; fa=endfull-blend;
          voi(ind)=(fb-ffff)/(fb-fa).*voivoiC + (ffff-fa)/(fb-fa).*voivoi;
        elseif (kk==5)   %%%%%right side
          %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
          %%%%%%%%%%%%%               fa                fb
          fa=endfull; fb=endfull+blend;          
          voi(ind)=(fb-ffff)/(fb-fa).*voivoi + (ffff-fa)/(fb-fa).*voivoiC;
          end
        end

    elseif ((LVF == 'B') | (LVF == 'b'))  
      if ((kk==1)|(kk==2))          %simply do cousin
        ymix=zeros(size(freqr)); 
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        %note we use k/kfar=1.0 instead of theratio(kk)
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,w_selfr,...
                           strenrt,stuff,layeramt,'V',tempbirn,1.0);
        voi(ind)=voivoi;

      elseif (kk == 3)               %do the fullmix
        %do full mixing at high resolution
        voivoi=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb);
        voif=boxint2(voivoi,boxcar);  

        %do the first mix second 
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);   
        voi1=boxint2(voivoi,boxcar);  
        voivoi=voi1;

        voi(ind)=(ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

      elseif ((kk==4)|(kk==5))   %do mixing = lorentz*ratio at output res
        %do full mixing at high resolution
        voivoi=full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb); 
        voif=voivoi;

        %do the first mix second 
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);   
        voi1=voivoi;
        voivoi=voi1;

        %now do the full/firstorder blend
        voila=(ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

        %do cousin
        ymix=zeros(size(freqr)); 
        %note we use k/kfar=1.0 instead of theratio(kk)
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        voici=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                 w_forr,w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        %now do the mixing/cousin blend!!!!!!!!!!!!!!!!
        %these are the start and stop points of current wavevector
        s1=ffff(1);   s2=ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend]=blender15um(band,prb,freqr,kk);
        if (kk==4)       %%%%%left side
          %%%%%%%%%%%%%%   COUS       |     BLEND        |       FULL
          %%%%%%%%%%%%%               fa                fb
          fb=endfull; fa=endfull-blend;
          voi(ind)=(fb-ffff)/(fb-fa).*voici + (ffff-fa)/(fb-fa).*voila;
        elseif (kk==5)   %%%%%right side
          %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
          %%%%%%%%%%%%%               fa                fb
          fa=endfull; fb=endfull+blend;          
          voi(ind)=(fb-ffff)/(fb-fa).*voila + (ffff-fa)/(fb-fa).*voici;
          end
        end
  
    elseif (IO == '1')
      if ((kk == 1)|(kk==2))     %do cousin
        ymix=zeros(size(freqr));  
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        %note we use k/kfar=1.0 instead of theratio(kk)
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,tempbirn,1.0);
        voi(ind)=voivoi;

      elseif (kk == 3)         %do mixing at high resolution
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        scum=voigtmix2(freqr,ffff,ymix,jr,...  
           temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);   
        voivoi=boxint2(scum,boxcar);
        voi(ind)=voivoi; 

      elseif ((kk == 4)|(kk==5))     %do cousin/mixing blend
        if quiet > 0
          fprintf(1,'doing cousin + linemix blend  .. ratio = 1.0 \n');
          end
        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);   
        voila=voivoi;                    %note we do not bioxcar avg here!!!
        ymix=zeros(size(freqr));  
        %note we use k/kfar=1.0 instead of theratio(kk)
        tempbirn=birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          tempbirn = 'c';
          end
        voici=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                w_forr,w_selfr,strenrt,stuff,layeramt,LVF,tempbirn,1.0);
        %these are the start and stop points of current wavevector
        s1=ffff(1);   s2=ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend]=blender15um(band,prb,freqr,kk);
        if (kk==4)       %%%%%left side
          %%%%%%%%%%%%%%   COUS       |     BLEND        |       FULL
          %%%%%%%%%%%%%               fa                fb
          fb=endfull; fa=endfull-blend;
          voi(ind)=(fb-ffff)/(fb-fa).*voici + (ffff-fa)/(fb-fa).*voila;
        elseif (kk==5)   %%%%%right side
          %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
          %%%%%%%%%%%%%               fa                fb
          fa=endfull; fb=endfull+blend;          
          voi(ind)=(fb-ffff)/(fb-fa).*voila + (ffff-fa)/(fb-fa).*voici;
          end
        end

    elseif (IO == '0')           %then birn = 'c' or 'n'
      if ((birn == 'b')|(birn == 'B'))
        error('need birn parameter = c or n')
        end
      if ((kk == 1)|(kk==2))     %do lorentz at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,1.0);
        voi(ind)=voivoi;
      elseif (kk == 3)     %do lorentz at high  res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                 w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,1.0);
        scum=boxint2(voivoi,boxcar);  
        voi(ind)=scum;

      elseif ((kk == 4)|(kk==5))     %do lorentz at output res
        %zeroth order line mixing calculation using voigt/lor lineshape  
        ymix=zeros(size(freqr));  
        voivoi=voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn,1.0);
        voi(ind)=voivoi;
        end             

      end         %if LVF = f or l or b

%    plot(outwave,exp(-voi),'+'); pause

    end           %if hinum(ind) > 0

  end             %outer for ii=1:5 loop

%check to make sure none of the terms are negative ....
voi=removeneg(outwave,voi);

%now smooth the transition from full mixing to RATIO mixing
%smooth over about 0.1 cm-1
df=abs(outwave(11)-outwave(1))/10; 
df=floor(0.05/df);
%df=floor(1/df);
df=1;
[voi]=smooth_full_ratio(outwave,voi,df,outnum,outstart,outstop);

%check again to make sure none of the terms are negative ....
voi=removeneg(outwave,voi);

fc=freqr;
cc=cputime-cc;

%if ((prb == 'R')|(prb == 'r'))
%  voi=zeros(size(voi));
%  end
