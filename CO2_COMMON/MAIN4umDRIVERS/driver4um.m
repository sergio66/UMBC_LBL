function [f,voi,fc,theRATIO]=...
        driver4um(temperature,ptotal,pself,layeramt,...  %gas layer params
                  band,prb,LVF,IO,birn,...               %spectra computation
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers
                  path_length,elower,jall,...            %from initial call
                  jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff_orig)

%this is for the 2500cm-1  (4 um band)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% set dofudge  = -1 for best results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%

%this section of code is identical for all bands of the 4 um mixing code
global quiet
global betafudge docfudge jmaxfudge

%%%%%%need slightly different parameters for "P" branch of 2350 ????

%%%might need to manipulate stuff_orig by fudging away, so make a copy of it
stuff=stuff_orig;

dofudge = +1;      %%%%this is to fit the data, or only use OLD RAL results
dofudge = -1;      %%%%this is RAL results + JJOHNS blend

dofudge = -1;
if dofudge > 0
  betafudge = 1.0;
  docfudge = 1.0;
  end

if (dofudge > 0)
  %%%this is when doing the CO2 JJOHNS fits
  if ((band == 2350) | (band == 2320) | (band == 2310)|(band==2351))
    disp(' ******  warning : dofudge == +1 in driver4um.m!!!!!!! *******');
    fprintf(1,'fudging away strong P,R pipi,sigsig,deltdelt ... \n');
    stuff.beta_pure=stuff_orig.beta_pure * betafudge;
    stuff.beta_for=stuff_orig.beta_for * betafudge;
    stuff.beta=stuff_orig.beta * betafudge;
    stuff.duration = stuff_orig.duration * docfudge;

    for ii = 1:3
      stuff.beta_pureJJ(ii)=stuff_orig.beta_pure * betafudge;
      stuff.beta_forJJ(ii)=stuff_orig.beta_for * betafudge;
      stuff.betaJJ(ii)=stuff_orig.beta * betafudge;
      stuff.durationJJ(ii) = stuff_orig.duration * docfudge;
      end

    end
   fprintf(1,'betafudge = %8.6f docfudge = %9.5f \n',betafudge,docfudge);
   fprintf(1,'and jmaxfudge = %8.6f \n',jmaxfudge);
fprintf(1,'beta = %8.6f dc  = %9.5f \n',stuff.beta,stuff.duration);
fprintf(1,'beta1= %8.6f dc1 = %9.5f \n',stuff.betaJJ(1),stuff.durationJJ(1));
fprintf(1,'beta2= %8.6f dc2 = %9.5f \n',stuff.betaJJ(2),stuff.durationJJ(2));
fprintf(1,'beta3= %8.6f dc3 = %9.5f \n',stuff.betaJJ(3),stuff.durationJJ(3));
  end

fprintf(1,'beta orig = %8.6f d of c = %9.5f \n',stuff.beta,stuff.duration);

cc = cputime;

f_hi = hiResFreq;
f_low = outwave;
f = outwave;

stuff.temperature = temperature;
stuff.prb         = prb;

[W_co2for,W_co2self] = wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuff);
[trans_ampl,population_t] = ...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuff);

W_plus = (pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

beta = stuff.beta;

W_plus = W_plus.*(beta+(1-beta)*eye(length(jr)));

if quiet > 0
  fprintf(1,'beta =  %8.6f d of c = %9.5f \n',beta,stuff.duration);
  end

%%%%%%%%%%%%%%%%%% adjust for weight effects!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%
iWeight = -1;         %best results!
iWeight = +1;         %not that much better, may even be worse

iWeight = -1;         %best results!

if iWeight > 0
  mass0 = 44;      %%%CO2 main isotope C-12 O-16 O-16
  mass1 = 44;
  end

%%%%%%%%%%%%%%%%%% end adjust for weight effects!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%

%%check doc values
if (stuff.duration < 1.0e-3)
  stuff.duration = 1.01e-3;
elseif (stuff.duration > 28.0e-3)
  stuff.duration = 27.0e-3;
  end
for ii = 1 : 3
  if (stuff.durationJJ(ii) < 1.0e-3)
    stuff.durationJJ(ii) = 1.01e-3;
  elseif (stuff.durationJJ(ii) > 28.0e-3)
    stuff.durationJJ(ii) = 27.0e-3;
    end
  end

%%check beta_pure values
if (stuff.beta_pure < 0)
  stuff.beta_pure = 0.001;
elseif (stuff.beta_pure > 1)
  stuff.beta_pure = 0.9999;
  end
for ii = 1 : 3
  %%check beta values
  if (stuff.beta_pureJJ(ii) < 0)
    stuff.beta_pureJJ(ii) = 0.001;
  elseif (stuff.beta_pureJJ(ii) > 1)
    stuff.beta_pureJJ(ii) = 0.9999;
    end
  end
 
%%check beta_for values
if (stuff.beta_for < 0)
  stuff.beta_for = 0.001;
elseif (stuff.beta_for > 1)
  stuff.beta_for = 0.9999;
  end
for ii = 1 : 3
  %%check beta values
  if (stuff.beta_forJJ(ii) < 0)
    stuff.beta_forJJ(ii) = 0.001;
  elseif (stuff.beta_forJJ(ii) > 1)
    stuff.beta_forJJ(ii) = 0.9999;
    end
  end

%%check beta values
if (stuff.beta < 0)
  stuff.beta = 0.001;
elseif (stuff.beta > 1)
  stuff.beta = 0.9999;
  end
for ii = 1 : 3
  %%check beta values
  if (stuff.betaJJ(ii) < 0)
    stuff.betaJJ(ii) = 0.001;
  elseif (stuff.betaJJ(ii) > 1)
    stuff.betaJJ(ii) = 0.9999;
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%if ((band == 2350) | (band == 2320) | (band == 2310))
if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
  fprintf(1,'beta = %8.6f d of c = %9.5f \n',stuff.beta,stuff.duration);

  stuffJJ1 = stuff;
  stuffJJ1.beta_pure       = stuff.beta_pureJJ(1);
  stuffJJ1.beta_for        = stuff.beta_forJJ(1);
  stuffJJ1.frequency_shift = stuff.frequency_shiftJJ(1);
  stuffJJ1.beta            = stuff.betaJJ(1);
  stuffJJ1.duration        = stuff.durationJJ(1);
  [W_co2forJJ1,W_co2selfJJ1]=wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuffJJ1);
  [trans_amplJJ1,population_tJJ1]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuffJJ1);
  W_plusJJ1 = ...
     (pself*W_co2selfJJ1+(ptotal-pself)*W_co2forJJ1)/stuffJJ1.pressure_ref;
  betaJJ1=stuffJJ1.beta;
  W_plusJJ1=W_plusJJ1.*(betaJJ1+(1-betaJJ1)*eye(length(jr)));
  fprintf(1,'beta-1= %8.6f dc-1= %9.5f\n',stuff.betaJJ(1),stuff.durationJJ(1));

  stuffJJ2 = stuff;
  stuffJJ2.beta_pure       = stuff.beta_pureJJ(2);
  stuffJJ2.beta_for        = stuff.beta_forJJ(2);
  stuffJJ2.frequency_shift = stuff.frequency_shiftJJ(2);
  stuffJJ2.beta            = stuff.betaJJ(2);
  stuffJJ2.duration        = stuff.durationJJ(2);
  [W_co2forJJ2,W_co2selfJJ2] = wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuffJJ2);
  [trans_amplJJ2,population_tJJ2]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuffJJ2);
  W_plusJJ2 = ...
     (pself*W_co2selfJJ2+(ptotal-pself)*W_co2forJJ2)/stuffJJ2.pressure_ref;
  betaJJ2=stuffJJ2.beta;
  W_plusJJ2=W_plusJJ2.*(betaJJ2+(1-betaJJ2)*eye(length(jr)));
  fprintf(1,'beta-2= %8.6f dc-2= %9.5f\n',stuff.betaJJ(2),stuff.durationJJ(2));

  stuffJJ3 = stuff;
  stuffJJ3.beta_pure       = stuff.beta_pureJJ(3);
  stuffJJ3.beta_for        = stuff.beta_forJJ(3);
  stuffJJ3.frequency_shift = stuff.frequency_shiftJJ(3);
  stuffJJ3.beta            = stuff.betaJJ(3);
  stuffJJ3.duration        = stuff.durationJJ(3);
  [W_co2forJJ3,W_co2selfJJ3] = wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuffJJ3);
  [trans_amplJJ3,population_tJJ3]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuffJJ3);
  W_plusJJ3 = ...
     (pself*W_co2selfJJ3+(ptotal-pself)*W_co2forJJ3)/stuffJJ3.pressure_ref;
  betaJJ3=stuffJJ3.beta;
  W_plusJJ3=W_plusJJ3.*(betaJJ3+(1-betaJJ3)*eye(length(jr)));
  fprintf(1,'beta-3= %8.6f dc-3= %9.5f\n',stuff.betaJJ(3),stuff.durationJJ(3));

elseif (band == 2351)
  fprintf(1,'beta = %8.6f d of c = %9.5f \n',stuff.beta,stuff.duration);

  stuffJJ1 = stuff;
  stuffJJ1.beta_pure       = stuff.beta_pureJJ(1);
  stuffJJ1.beta_for        = stuff.beta_forJJ(1);
  stuffJJ1.frequency_shift = stuff.frequency_shiftJJ(1);
  stuffJJ1.beta            = stuff.betaJJ(1);
  stuffJJ1.duration        = stuff.durationJJ(1);
  [W_co2forJJ1,W_co2selfJJ1] = wfunco2er(jr,elower,elowerr,w_selfr,w_forr,...
                                   band,jall,temperature,stuffJJ1);
  [trans_amplJJ1,population_tJJ1]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuffJJ1);
  W_plusJJ1 = ...
     (pself*W_co2selfJJ1+(ptotal-pself)*W_co2forJJ1)/stuffJJ1.pressure_ref;
  betaJJ1=stuffJJ1.beta;
  W_plusJJ1=W_plusJJ1.*(betaJJ1+(1-betaJJ1)*eye(length(jr)));
  fprintf(1,'beta-1= %8.6f dc-1= %9.5f\n',stuff.betaJJ(1),stuff.durationJJ(1));
  end

ratio=doratio(population_t,trans_ampl,W_plus);
if (ratio <= 0.0) 
  ratio=1.0e-10;         %if the sum rule give -ve numbers
  end
theRATIO=ratio;

%now find out where the boundaries to do full mixing and just doing 
%ratio to Lorentz, ie
[hinum,histart,histop,outnum,outstart,outstop,theratio]=...
   tempratio4um(band,temperature,boxcar,f_hi,f_low,prb,freqr,ratio,...
                       strenrt,w_forr,w_selfr,jr,stuff);

voi=zeros(size(outwave));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk=1:5
  theratio(kk)=ratio;

  if hinum(kk) > 0
    ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
    %%if ((band == 2350) | (band == 2320) | (band == 2310))
    if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
      ymixJJ1 = ...
        y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,stuffJJ1,betaJJ1);
      ymixJJ2 = ...
        y1s(jr,freqr,elowerr,strenr,W_plusJJ2,trans_amplJJ2,stuffJJ2,betaJJ2);
      ymixJJ3 = ...
        y1s(jr,freqr,elowerr,strenr,W_plusJJ3,trans_amplJJ3,stuffJJ3,betaJJ3);
    elseif (band == 2351)
      ymixJJ1 = ...
        y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,stuffJJ1,betaJJ1);
      end

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

    if quiet > 0
      fprintf(1,'--> kk= %3i freqs = %10.5f %10.5f \n',kk,min(ffff),max(ffff));
      end

    if ((LVF == 'F') | (LVF == 'f'))             %XXXXXXXXXXXXXX
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
        ymix = zeros(size(freqr)); 
        voivoi = voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,...
            w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        voi(ind) = voivoi;
        %plot(ffff,exp(-voivoi)); title('kk = 1,2'); pause(0.1);

      elseif (kk == 3)          %do full mixing at high resolution
        scumRAL = full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                   trans_ampl,population_t,stuff,birn,theratio(kk),...
                   strenrt,layeramt,ymix,band,prb);
        %%if ((band == 2350)  | (band == 2320) | (band == 2310))
        if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
          scumJJ1 = ...
            full2(freqr,ffff,W_plusJJ1,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ1,population_tJJ1,stuffJJ1,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ1,band,prb);
          scumJJ2 = ...
            full2(freqr,ffff,W_plusJJ2,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ2,population_tJJ2,stuffJJ2,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ2,band,prb);
          scumJJ3 = ...
            full2(freqr,ffff,W_plusJJ3,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ3,population_tJJ3,stuffJJ3,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ3,band,prb);
          scum = ...
            ral_jjohn_blend3(band,prb,ffff,scumRAL,scumJJ1,scumJJ2,scumJJ3);
        elseif (band == 2351)
          scumJJ1=full2(freqr,ffff,W_plusJJ1,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ1,population_tJJ1,stuffJJ1,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ1,band,prb);
          scum = ral_jjohn_blend(band,prb,freqr,strenrt,ffff,scumRAL,scumJJ1);
        else
          scum = scumRAL;
          end
        voivoi = boxint2(scum,boxcar);  
        %plot(ffff,exp(-scum)); title('LVF = f, kk=3 == high res'); pause(0.1)
        voi(ind) = voivoi;

      elseif ((kk==4)|(kk==5))      %do cousin + linemixing blend
        if quiet > 0
          fprintf(1,'doing cousin + linemix blend .. ratio = 1.0 \n');
          end
        voivoi = full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb); 
        ymix = zeros(size(freqr)); 
        %note we use k/kfar = 1.0 instead of theratio(kk)
        tempbirn = birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          tempbirn = 'c';
          end
        voivoiC = voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,...
            w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        %these are the start and stop points of current wavevector
        s1 = ffff(1);   s2 = ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend] = blender4um(band,prb,freqr,jr,kk);
        if (kk==4)       %%%%%left side
          %%%%%%%%%%%%%%   COUS       |     BLEND        |       FULL
          %%%%%%%%%%%%%               fa                fb
          fb = endfull; fa = endfull-blend;
          voi(ind) = (fb-ffff)/(fb-fa).*voivoiC + (ffff-fa)/(fb-fa).*voivoi;
        elseif (kk==5)   %%%%%right side
          %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
          %%%%%%%%%%%%%               fa                fb
          fa = endfull; fb = endfull+blend;          
          voi(ind) = (fb-ffff)/(fb-fa).*voivoi + (ffff-fa)/(fb-fa).*voivoiC;
          end
        %plot(ffff,exp(-voivoi)); title('kk = 4,5'); pause(0.1);
        end

    elseif ((LVF == 'B') | (LVF == 'b'))             %XXXXXXXXXXXXXX
      if ((kk==1)|(kk==2))          %simply do cousin
        ymix = zeros(size(freqr)); 
        tempbirn = birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        %note we use k/kfar=1.0 instead of theratio(kk)
       voivoi = voigtmixRatio(freqr,ffff,ymix,jr,temperature,w_forr,w_selfr,...
                           strenrt,stuff,layeramt,'V',tempbirn,1.0);
        voi(ind) = voivoi;

      elseif (kk == 3)               %do the fullmix
        %do full mixing at high resolution
        scumRAL = full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb);
        %if ((band == 2350)  | (band == 2320) | (band == 2310))
        if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
          scumJJ1 = ...
            full2(freqr,ffff,W_plusJJ1,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ1,population_tJJ1,stuffJJ1,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ1,band,prb);
          scumJJ2 = ...
            full2(freqr,ffff,W_plusJJ2,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ2,population_tJJ2,stuffJJ2,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ2,band,prb);
          scumJJ3 = ...
            full2(freqr,ffff,W_plusJJ3,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ3,population_tJJ3,stuffJJ3,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ3,band,prb);
          scum = ...
            ral_jjohn_blend3(band,prb,ffff,scumRAL,scumJJ1,scumJJ2,scumJJ3);
        elseif (band == 2351) 
         scumJJ1 = full2(freqr,ffff,W_plusJJ1,jr,w_selfr,w_forr,temperature,...
                   trans_amplJJ1,population_tJJ1,stuffJJ1,birn,theratio(kk),...
                   strenrt,layeramt,ymixJJ1,band,prb);
         scum  =  ral_jjohn_blend(band,prb,freqr,strenrt,ffff,scumRAL,scumJJ1);
        else
          scum = scumRAL;
          end
        %plot(ffff,scum); title('LVF = b so fullmix kk=3 == high res'); pause
        voif=boxint2(scum,boxcar);  

%%        %do the first mix second 
%%        ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
%%        voivoi=voigtmix2(freqr,ffff,ymix,jr,...  
%%              temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);
%%      %plot(ffff,scum); title('LVF = b so first mix kk=3 == high res'); pause
%%        voi1=boxint2(voivoi,boxcar);  
%%        voivoi=voi1;
        ymixRAL=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        scumRAL=voigtmix2(freqr,ffff,ymixRAL,jr,...  
           temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);   
        %if ((band == 2350)  | (band == 2320) | (band == 2310))
        if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
          ymixJJ1 = y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,...
                      stuffJJ1,betaJJ1);
          scumJJ1 = voigtmix2(freqr,ffff,ymixJJ1,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ1,layeramt,LVF,birn);   
          ymixJJ2 = y1s(jr,freqr,elowerr,strenr,W_plusJJ2,trans_amplJJ2,...
                      stuffJJ2,betaJJ2);
          scumJJ2 = voigtmix2(freqr,ffff,ymixJJ2,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ2,layeramt,LVF,birn);   
          ymixJJ3 = y1s(jr,freqr,elowerr,strenr,W_plusJJ3,trans_amplJJ3,...
                      stuffJJ3,betaJJ3);
          scumJJ3 = voigtmix2(freqr,ffff,ymixJJ3,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ3,layeramt,LVF,birn);   

        scum = ral_jjohn_blend3(band,prb,ffff,scumRAL,scumJJ1,scumJJ2,scumJJ3);
        elseif (band == 2351)
          ymixJJ1 = y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,...
                      stuffJJ1,betaJJ1);
          scumJJ1 = voigtmix2(freqr,ffff,ymixJJ1,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ1,layeramt,LVF,birn);   
          scum = ral_jjohn_blend(band,prb,freqr,strenrt,ffff,scumRAL,scumJJ1);
        else
          scum = scumRAL;
          end
        %plot(ffff,scum); title('IO=1,blend kk=3 == high res'); pause
        voi1=boxint2(scum,boxcar);
        voivoi = voi1; 

        voi(ind) = (ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

      elseif ((kk==4)|(kk==5))   %do mixing = lorentz*ratio at output res
        %do full mixing at high resolution
        voivoi = full2(freqr,ffff,W_plus,jr,w_selfr,w_forr,temperature,...
                     trans_ampl,population_t,stuff,birn,theratio(kk),...
                     strenrt,layeramt,ymix,band,prb);
        voif = voivoi;
        %do the first mix second 
        ymix = y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        voivoi = voigtmix2(freqr,ffff,ymix,jr,...  
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);   
        voi1 = voivoi;
        voivoi = voi1;

        %now do the full/firstorder blend
        voila = (ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 

        %do cousin
        ymix = zeros(size(freqr)); 
        %note we use k/kfar = 1.0 instead of theratio(kk)
        tempbirn = birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        voici = voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
                 w_forr,w_selfr,strenrt,stuff,layeramt,'V',tempbirn,1.0);
        %these are the start and stop points of current wavevector
        s1 = ffff(1);   s2 = ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend] = blender4um(band,prb,freqr,jr,kk);
        if (kk==4)       %%%%%left side
          %%%%%%%%%%%%%%   COUS       |     BLEND        |       FULL
          %%%%%%%%%%%%%               fa                fb
          fb = endfull; fa = endfull-blend;
          voi(ind) = (fb-ffff)/(fb-fa).*voici + (ffff-fa)/(fb-fa).*voila;
        elseif (kk==5)   %%%%%right side
          %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
          %%%%%%%%%%%%%               fa                fb
          fa = endfull; fb = endfull+blend;          
          voi(ind) = (fb-ffff)/(fb-fa).*voila + (ffff-fa)/(fb-fa).*voici;
          end
        end
  
    elseif (IO == '1')                  %XXXXXXXXXXXXXX
      if ((kk == 1)|(kk==2))     %do cousin
        ymix = zeros(size(freqr));  
        tempbirn = birn;
        if (tempbirn == 'b')|(tempbirn == 'B')
          if quiet > 0
            fprintf(1,'doing cousin instead of lorentz .. ratio = 1.0 \n');
            end
          tempbirn = 'c';
          end
        %note we use k/kfar=1.0 instead of theratio(kk)
        voivoi = voigtmixRatio(freqr,ffff,ymix,jr,temperature,...
              w_forr,w_selfr,strenrt,stuff,layeramt,LVF,tempbirn,1.0);
        voi(ind) = voivoi;

      elseif (kk == 3)         %do mixing at high resolution
        ymixRAL = y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
        scumRAL = voigtmix2(freqr,ffff,ymix,jr,...  
           temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);   

        %if ((band  == 2350)  | (band == 2320) | (band == 2310))
        if ((band == 2350) & ((prb == 'R') | (prb == 'r')))
          ymixJJ1=y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,...
                      stuffJJ1,betaJJ1);
          scumJJ1=voigtmix2(freqr,ffff,ymixJJ1,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ1,layeramt,LVF,birn);   
          ymixJJ2=y1s(jr,freqr,elowerr,strenr,W_plusJJ2,trans_amplJJ2,...
                      stuffJJ2,betaJJ2);
          scumJJ2=voigtmix2(freqr,ffff,ymixJJ2,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ2,layeramt,LVF,birn);   
          ymixJJ3=y1s(jr,freqr,elowerr,strenr,W_plusJJ3,trans_amplJJ3,...
                      stuffJJ3,betaJJ3);
          scumJJ3=voigtmix2(freqr,ffff,ymixJJ3,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ3,layeramt,LVF,birn);   

          scum=ral_jjohn_blend3(band,prb,ffff,scumRAL,scumJJ1,scumJJ2,scumJJ3);
        elseif (band == 2351)
          ymixJJ1=y1s(jr,freqr,elowerr,strenr,W_plusJJ1,trans_amplJJ1,...
                      stuffJJ1,betaJJ1);
          scumJJ1=voigtmix2(freqr,ffff,ymixJJ1,jr,...  
             temperature,w_forr,w_selfr,strenrt,stuffJJ1,layeramt,LVF,birn);   
          scum = ral_jjohn_blend(band,prb,freqr,strenrt,ffff,scumRAL,scumJJ1);
        else
          scum = scumRAL;
          end
        %plot(ffff,scum); title('IO=1,no blend kk=3 == high res'); pause
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
        [endfull,blend]=blender4um(band,prb,freqr,jr,kk);
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

    elseif (IO == '0')           %then birn = 'c' or 'n'  %XXXXXXXXXXXXXX
      fprintf(1,'hmmm either doing cousin or nothing! %3i %s \n',kk,birn)
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
        %plot(ffff,scum); title('LVF = 0 so no mix kk=3 == high res'); pause
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

    if quiet > 0    
      plot(outwave,exp(-voi)); pause(0.1);
      end

    end           %if hinum(ind) > 0


  end             %outer for ii=1:5 loop

%check to make sure none of the terms are negative ....

voi=removeneg(outwave,voi);

%now smooth the transition from full mixing to RATIO mixing
%smooth over about 0.1 cm-1
df=abs(outwave(11)-outwave(1))/10; 
%df=floor(0.05/df);
df=floor(0.005/df);
df=1;
df=0;
[voi]=smooth_full_ratio(outwave,voi,df,outnum,outstart,outstop);

%check again to make sure none of the terms are negative ....
voi=removeneg(outwave,voi);

fc=freqr;
cc=cputime-cc;

