function [f,voi,fc,theRATIO]=...
   driver15umPR_SIMPLE(temperature,ptotal,pself,layeramt,...  %gas layer params
                  band,prb,LVF,IO,birn,...               %spectra computation
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers
                  path_length,...            %from initial call
                  jq,w_forq,w_selfq,freqq,strenqt,stuff)

%this is the file used for computations if k=klor*constant

%this is for the 700cm-1  (15 um band)
%
%  ----------------|-----------------------------------------|------------- 
%                  f1                                        f2 
%   k=klor*ratio                 k=kfull                        k=klor*ratio 
% 
%         I                        III                            II 
% 
 
if ((prb == 'q') | (prb == 'Q'))
  prb
  error('this driver is for the PR sigpi,deltpi mixing calcs!!')
  end
  
%if (abs(sigORdelt) ~=1) 
%  sigORdelt
%  error('this driver is for the PR sigpi,deltpi mixing calcs!!')
%  end
  

%this section of code is identical for Q bands of the 15 um mixing code
cc=cputime;

f_hi=hiResFreq;
f_low=outwave;
f=outwave;

%%% other bands have "simple" linemixing
%now find out where the boundaries to do full mixing and just doing  
%ratio to Lorentz, ie 

if quiet > 0
  fprintf(1,'directly computing lineshape for PR band = %3i\n',band);
  fprintf(1,'if cousin not asked for, uses k/klor = %8.6f \n',theRATIO);
  end

ratio=theRATIO;
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
      
    voivoi=voigt_mix_simple2(freqq,ffff,ymix,jq,temperature,w_forq,...
         w_selfq,strenqt,stuff,layeramt,LVFtemp,IOtemp,birn,theRATIO);   
    if (kk == 3)
      scum=boxint2(voivoi,boxcar);
      voi(ind)=scum;
    else
      voi(ind)=voivoi;
      end

    end             %if hinum(kk) > 0
  end               %for kk=1,3

fc=freqq;
cc=cputime-cc;
