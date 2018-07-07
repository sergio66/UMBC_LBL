function [f,voi,fc,theRATIO]=...
   driver15umPR(temperature,ptotal,pself,layeramt,...  %gas layer params
                  band,prb,LVF,IO,birn,...               %spectra computation
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers
                  path_length,...            %from initial call
                  jq,w_forq,w_selfq,freqq,strenqt,stuff)

%this is the file used for computations if k=klor*constant

%this is for the 15 um PR mixing
% 
%  -------|------------|----------------------------|-------------|---------- 
%        f1-1          f1                          f2           f2+1 
%   
%  cousin   full/cousin          full               full/cousin     cousin     
% 
%   I          IV                 III                   V            II 
% 

global quiet
 
if ((prb == 'q') | (prb == 'Q'))
  prb
  error('this driver is for the 15um PR mixing calcs!!')
  end

%this section of code is identical for Q bands of the 15 um mixing code
cc=cputime;

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

%%% other bands have "simple" linemixing
%now find out where the boundaries to do full mixing and just doing  
%ratio to Lorentz, ie 

if quiet > 0
  fprintf(1,'directly computing lineshape for 15 um PR band = %3i\n',band);
  fprintf(1,'if cousin not asked for, uses k/klor = %8.6f \n',theRATIO);
  end
ratio=theRATIO;
fprintf(1,' PR mixing : using ratio = %8.6f for %s branch \n',theRATIO,prb);

stuff.temperature = temperature;
stuff.prb         = prb;

[hinum,histart,histop,outnum,outstart,outstop,theratio] = ...
     tempratio15um(band,temperature,boxcar,f_hi,f_low,prb,freqq,ratio,...
                   strenqt,w_forq,w_selfq,jq,stuff);

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

for kk=1:5
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

    if (kk == 3)      %do the computations with simple line mix
      voivoi=voigt_mix_simple2(freqq,ffff,ymix,jq,temperature,w_forq,...
         w_selfq,strenqt,stuff,layeramt,LVFtemp,IOtemp,birn,theRATIO);   
      scum=boxint2(voivoi,boxcar);
      voi(ind)=scum;
      end

    if ((kk == 1)|(kk==2))      %do the computations with COUSIN only
      voivoiC=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
         w_selfq,strenqt,stuff,layeramt,'V','c',1.0);
      voi(ind)=voivoiC;
      end

    if ((kk == 4)|(kk==5))   %do the computations, blending COUSIN and LINE MIX
      voivoiC=voigtmixRatio(freqq,ffff,ymix,jq,temperature,w_forq,...
         w_selfq,strenqt,stuff,layeramt,'V','c',1.0);   
      voivoi =voigt_mix_simple2(freqq,ffff,ymix,jq,temperature,w_forq,...
         w_selfq,strenqt,stuff,layeramt,LVFtemp,IOtemp,birn,theRATIO);   
      iOld = +1;     %%this is old, do not use
      iOld = -1;     %%this was the code in Summer 2007
      if iOld > 0
        %%%%old code
        slpe=(ffff(length(ffff))-ffff(1)); 
        s1=ffff(1);   s2=ffff(length(ffff)); 
        if (band < s1)
          voi(ind)=((voivoiC.*(ffff-s1) + voivoi.*(s2-ffff)))/slpe; 
        else
          voi(ind)=((voivoi.*(ffff-s1) + voivoiC.*(s2-ffff)))/slpe; 
          end
      else
        %these are the start and stop points of current wavevector
        s1=ffff(1);   s2=ffff(length(ffff));
        %endfull tells where the end of full mixing is, and blend is the
        %width of the blending region
        [endfull,blend]=blender15um(band,prb,freqq,kk);
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
      end
    end             %if hinum(kk) > 0
  end               %for kk=1,5

%check to make sure none of the terms are negative ....
voi=removeneg(outwave,voi);

%now smooth the transition from full mixing to RATIO mixing
%smooth over about 0.1 cm-1
df=abs(outwave(11)-outwave(1))/10; 
df=floor(0.05/df);
%df=floor(1/df);
[voi]=smooth_full_ratio(outwave,voi,df,outnum,outstart,outstop);

%check again to make sure none of the terms are negative ....
voi=removeneg(outwave,voi);

fc=freqq;
cc=cputime-cc;
