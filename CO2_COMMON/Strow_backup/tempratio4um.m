function [hinum,histart,histop,outnum,outstart,outstop,theratio]=... 
          tempratio4um(band,temp,boxcar,hifreq,outwave,prb,freqr,INratio,...
                       strenrt,w_forr,w_selfr,jr,stuff);

%%%%%%% WARNING : IF YOU MAKE CHANGES HERE, MAKE CHANGES IN BLENDER4UM.M
%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%

% INratio is what comes from Strow's paper ... it could be quite small
% for the 4 um band, so we have to go further away before we can safely
% say      k=klor*ratio

%this function tells you how the full/lorentz ratio varies with temp 
%the answers are 1d matrices of 3 elements each :  
%    kk=1 : f < f1   -------- do cousin
%    kk=2 : f > f2   -------- do cousin
%    kk=3 : (f > f1)&(f < f2)         -------- do full mixing 
%    kk=4 : f1-1 < f < f1  -------- do cousin/full mixing blend
%    kk=5 : f2 < f < f2+1  -------- do cousin/full mixing blend
%

% for the HIGH resolution wavevector f (eg 0.0005 wavenumber spacing) : since 
% there is a boxcar average done, this number is used in the computation of 
% the start/stop points 
%   hinum tells you how many elements are in <hifreq> region kk 
%   histart tells you where index starts for region of interest 
%   histop  tells you where index stops  for region of interest 

% for the lower OUTPUT resolution wavevector <outwave> 
% (eg 0.0025 wavenumber spacing) 
%   outnum tells you how many elements are in freq region kk 
%   outstart tells you where index starts for region of interest 
%   outstop  tells you where index stops  for region of interest 

% theratio gives you the full/lorentz ratio for BOTH resolutions 
%the data below is in the form [wavenumber ratio(1) .. ratio(7)] where 
% theratio = full/lor for following 7 temperatures temp(1):temp(7) =  ... 
%                                    150 200 .. 400 
% actually, theratio is usually set to +/-1.0, as it is accurately computed
% using Strows eqn

% hifreq  = input wavenumber (at the high resolution (before boxcar avg))

%%%%%%%%%%%%%%%%%%%%%% the useless stuff !!!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% thus, the bottom line is that output parameter theratio is set to +/- 1  %%
%% in what follows, only matr(2,1) and matr(3,1) are used, as               %%
%% they set the boundaries iupto which full/full+cousin/cousin is done      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% possible problems %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if the input wavevector "hifreq" is too coarse, then the code might find %%
%% indices that are not commensurate with boxcar everaging. eg if invector  %%
%% starts at 2150, we gave P 2350 branch, whose min line is at 2230, then   %%
%% matr(3,1)=min(freqr)-80 == 2150!!! so outstop(1)  == 1 ... so there will %%
%% be problems when we try to find histop, as we have to set the "lolo" by  %%
%% going (boxcar-1) points to the left of this value !!!!!!                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global quiet p2311_21_jmax pr2350_jmax pr2351_jmax

epep=abs(hifreq(length(hifreq))-hifreq(1))/(length(hifreq)-1)/10.0;
bbx=(boxcar-1)/2;
 
%set up dummy variables if necessary
if (nargin==5) 
  prb='Q'; freqr=[1]; 
elseif (nargin==6) 
  freqr=[1]; 
  end 

find_regions; 

%%%%%%% do the lower OUTPUT resolution stuff first %%%%%%%%%%%%%%%%%%%%%%
%now find the indices for 5 regions of interest
%this is for the OUTPUT resolution

if (band ~= 2350) 
  iiLessOut=find(outwave < matr(2,1));                                  %%Reg 1
    numLessOut=length(iiLessOut);
  iiMoreOut=find(outwave >= matr(3,1));                                 %%Reg 2
   numMoreOut=length(iiMoreOut);
  iiNumOut =find((outwave >= matr(2,1)+1) & (outwave < matr(3,1)-1));   %%Reg 3
    numOut=length(iiNumOut);
  iiLessMidOut =find((outwave >= matr(2,1)) & (outwave < matr(2,1)+1)); %%Reg 4
    numLessMidOut=length(iiLessMidOut); 
  iiMoreMidOut =find((outwave >= matr(3,1)-1) & (outwave < matr(3,1))); %%Reg 5
    numMoreMidOut=length(iiMoreMidOut);
else  %%%%%make the blending region larger for the strong 2350 band !!!!!!!!!!!
  iiLessOut=find(outwave < matr(2,1));                                  %%Reg 1
    numLessOut=length(iiLessOut);
  iiMoreOut=find(outwave >= matr(3,1));                                 %%Reg 2
   numMoreOut=length(iiMoreOut);
  iiNumOut =find((outwave >= matr(2,1)+20) & (outwave < matr(3,1)-20)); %%Reg 3
    numOut=length(iiNumOut);
  iiLessMidOut =find((outwave >= matr(2,1)) & (outwave < matr(2,1)+20));%%Reg 4
    numLessMidOut=length(iiLessMidOut); 
  iiMoreMidOut =find((outwave >= matr(3,1)-20) & (outwave < matr(3,1)));%%Reg 5
    numMoreMidOut=length(iiMoreMidOut);
  end

for kk = 1:5
  %find the start and stop freqs of the regions of interest
  %this is for lower res outwavevector
  if (kk == 1)
    indexOut=iiLessOut;
    tempfreqOut=outwave(iiLessOut);  num=numLessOut;
  elseif (kk == 2)
    indexOut=iiMoreOut;
    tempfreqOut=outwave(iiMoreOut);  num=numMoreOut;
  elseif (kk == 3)
    indexOut=iiNumOut;
    tempfreqOut=outwave(iiNumOut);   num=numOut;
  elseif (kk == 4)
    indexOut=iiLessMidOut;
    tempfreqOut=outwave(iiLessMidOut);   num=numLessMidOut;
  elseif (kk == 5)
    indexOut=iiMoreMidOut;
    tempfreqOut=outwave(iiMoreMidOut);   num=numMoreMidOut;
    end

  if (num > 0)
    frstartOut=tempfreqOut(1);  frendOut=tempfreqOut(length(tempfreqOut));
    ratio=1.0;
    outnum(kk)=num;
    outratio(kk)=ratio;
    if (kk == 3)      %this is the region of full mixing
      outratio(kk)=1.0;
      end
    outstart(kk)=indexOut(1);
    outstop(kk)=indexOut(length(indexOut));

  else
    outnum(kk)=0;
    outratio(kk)=-1.0;
    outstart(kk)=-1;
    outstop(kk)=-1;
    end  %if loop

  end  %for loop

%%%%%%%%% do the HIGH resolution stuff second %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now find the indices for 5 regions of interest
%this is for the HIGH resolution
if (band ~= 2350) 
  iiLessHi=find(hifreq < matr(2,1));                                %%Reg 1
    numLessHi=length(iiLessHi);
  iiMoreHi=find(hifreq >= matr(3,1));                               %%Reg2
    numMoreHi=length(iiMoreHi);
  iiNumHi =find((hifreq >= matr(2,1)+1) & (hifreq < matr(3,1)-1));  %%Reg 3
    numHi=length(iiNumHi);
  iiLessMidHi =find((hifreq >= matr(2,1)) & (hifreq < matr(2,1)+1));%%Reg 4
    numLessMidHi=length(iiLessMidHi);
  iiMoreMidHi =find((hifreq >= matr(3,1)-1) & (hifreq < matr(3,1)));%%Reg 5
    numMoreMidHi=length(iiMoreMidHi);
else          %%%%%%wider blend regions!!!
  iiLessHi=find(hifreq < matr(2,1));                                %%Reg 1
    numLessHi=length(iiLessHi);   
  iiMoreHi=find(hifreq >= matr(3,1));                               %%Reg2
    numMoreHi=length(iiMoreHi);
  iiNumHi =find((hifreq >= matr(2,1)+20) & (hifreq < matr(3,1)-20));%%Reg 3
    numHi=length(iiNumHi);
  iiLessMidHi =find((hifreq >= matr(2,1)) & (hifreq < matr(2,1)+20));%%Reg 4
    numLessMidHi=length(iiLessMidHi);
  iiMoreMidHi =find((hifreq >= matr(3,1)-20) & (hifreq < matr(3,1)));%%Reg 5
    numMoreMidHi=length(iiMoreMidHi);
  end

for kk = 1:5
  %find the start and stop freqs of the regions of interest
  %this is for high res wavevector
  if (kk == 1)
    indexHi=iiLessHi;
    tempfreqHi=hifreq(iiLessHi);  num=numLessHi;
  elseif (kk == 2)
    indexHi=iiMoreHi;
    tempfreqHi=hifreq(iiMoreHi);  num=numMoreHi;
  elseif (kk == 3)
    indexHi=iiNumHi;
    tempfreqHi=hifreq(iiNumHi);    num=numHi;
  elseif (kk == 4)
    indexHi=iiLessMidHi;
    tempfreqHi=hifreq(iiLessMidHi);   num=numLessMidHi;
  elseif (kk == 5)
    indexHi=iiMoreMidHi;
    tempfreqHi=hifreq(iiMoreMidHi);   num=numMoreMidHi;
    end

  if (num > 0)
    frstartHi=tempfreqHi(1);  frendHi=tempfreqHi(length(tempfreqHi));
    ratio = 1.0;
    %this is ignoring the boxcar
    hinum(kk)=num;
    theratio(kk)=ratio;

    if (kk == 3)      %this is the region of full mixing
      theratio(kk)=1.0;
      end

    histart(kk)=indexHi(1);
    histop(kk)=indexHi(length(indexHi));

    %now add on effects of boxcar
    %do the START of the vector ... two lines below are commented because 
    %outstart(kk) < 0 is possible
    %outSTART=outwave(outstart(kk));
    %hihi=boxint2(hifreq(histart(kk):histart(kk)+boxcar-1),boxcar);

    if (outstart(kk) == 1)       %then we need histart(kk) == 1 also
      if (histart(kk) ~= 1)
        extra=histart(kk)-1;
        histart(kk) = 1;
        hinum(kk)=hinum(kk)+extra;
        end
      end
 
    if (outstart(kk) > 1)      %check to make sure things are ok
      outSTART=outwave(outstart(kk));
      hihi=boxint2(hifreq(histart(kk):histart(kk)+boxcar-1),boxcar);
      dada=hihi-outSTART;
      if (abs(dada) > epep)    %%%oh boy : have to adjust!
        truth= -1;
        histart(kk)=histart(kk)-boxcar;
        hinum(kk)=hinum(kk)+boxcar; %ADDED on more points to start of vector
        while (truth < 0)
          hihi=boxint2(hifreq(histart(kk):histart(kk)+boxcar-1),boxcar);
          dada=hihi-outSTART;            
          if (abs(dada) > epep)
            histart(kk)=histart(kk)+1;
            hinum(kk)=hinum(kk)-1;  %taking AWAY points from start of vector
          else
            truth = +1;
            end
          end
        end
      end

    if (outstart(kk) < 0)      %check to make sure things are ok
      histart(kk)=-1;
      histop(kk)=-1;
      hinum(kk)=0;       %NO POINTS in vector
      end

    %do the END of the vector ... two lines below are commented because 
    %outstop(kk) < 0 is possible
    %outSTOP=outwave(outstop(kk));
    %lolo=boxint2(hifreq(histop(kk)-(boxcar-1):histop(kk)),boxcar);

    if (outstop(kk) == length(outwave))
      %then we need histop(kk) == length(hifreq) also
      if (histop(kk) ~= length(hifreq))
        extra=length(hifreq)-histop(kk);
        histop(kk) = length(hifreq);
        hinum(kk)=hinum(kk)+extra;
        end
      end
 
    if ( (outstop(kk) < length(outwave)) & (outstop(kk) > 0))
      %check to make sure things are ok
      outSTOP=outwave(outstop(kk));
      lolo=boxint2(hifreq(histop(kk)-(boxcar-1):histop(kk)),boxcar);
      dada=lolo-outSTOP;
      if (abs(dada) > epep)    %%%oh boy : have to adjust!
        truth= -1;
        histop(kk)=histop(kk)-boxcar;
        hinum(kk)=hinum(kk)-boxcar;
        while (truth < 0)
          lolo=boxint2(hifreq(histop(kk)-(boxcar-1):histop(kk)),boxcar);
          dada=lolo-outSTOP;            
          if ((abs(dada) > epep) & (histop(kk) < length(hifreq)))
            histop(kk)=histop(kk)+1;
            hinum(kk)=hinum(kk)+1;
          elseif ((abs(dada) > epep) & (histop(kk) == length(hifreq)))
            error('cannot go beyond end of vector!!!');
          else
            truth = +1;
            end
          end
        end
      end

    if (outstop(kk) < 0)      %check to make sure things are ok
      histart(kk)=-1;
      histop(kk)=-1;
      hinum(kk)=0;       %NO POINTS in vector
      end


%this next set of if-thens is just for debug purposes!
    if (outstart(kk) > 0)
      outSTART=outwave(outstart(kk));  
    else 
       outSTART=-1.00;
       end
    if (histart(kk) > 0)
      hihi=boxint2(hifreq(histart(kk):histart(kk)+boxcar-1),boxcar);
    else 
       hihi=-1.00;
       end

    if (outstop(kk) > 0)
      outSTOP=outwave(outstop(kk));
    else 
       outSTOP=-1.00;
       end
    if (histop(kk) > 0)
      lolo=boxint2(hifreq(histop(kk)-(boxcar-1):histop(kk)),boxcar);
    else 
       lolo=-1.00;
       end

  else  %%%%%%if num <= 0
    hinum(kk)=0;
    theratio(kk)=-1.0;
    histart(kk)=-1;
    histop(kk)=-1;
    end  %if loop

  end  %for loop

%%%%%%%%%%%%%%%% print results to see if they make sense
quiet0=quiet;

%%%%%%%%%quiet=1;

pprr=quiet;
if pprr > 0
  fprintf(1,'printing out stuff in tempratio4um \n');

  fprintf(1,'matr(2,1) matr(3,1) = %10.6f %10.6f \n',matr(2,1),matr(3,1));

  summ=0;
  for ii=1:5     
    if (outnum(ii) > 0)
      doh=[ii  outstart(ii) outstop(ii) ...
           outwave(outstart(ii)) outwave(outstop(ii)) outnum(ii)];
      summ=summ+outnum(ii);
      fprintf(1,'%3i %6i %6i %9.7e  %9.7e  %6i \n',doh);
      end;
    end
  fprintf(1,'sum of pts, length(outwave) = %10i %10i \n',summ,length(outwave));

  summ=0;
  for ii=1:5     
    if (hinum(ii) > 0)
      doh=[ii   histart(ii) histop(ii) ...
           hifreq(histart(ii)) hifreq(histop(ii)) hinum(ii)];
      summ=summ+hinum(ii);
      fprintf(1,'%3i %6i %6i %9.7e  %9.7e  %6i \n',doh);
      end;
    end
  fprintf(1,'sum of pts, length(hiwave) = %10i %10i \n',summ,length(hifreq));

  end

if quiet > 0
  format short e
  fprintf(1,'freqr min and max   are %10.6f %10.6f \n',min(freqr),max(freqr))
  fprintf(1,'alloted min and max are %10.6f %10.6f \n',matr(2,1),matr(3,1))
  kk=1;
  if hinum(kk) > 0
    a1=hifreq(histart(kk));
    a2= hifreq(histop(kk));
    fprintf(1,'kk = %3i start, stop = %10.6f %10.6f \n',kk,a1,a2);
    end
  kk=4;
  if hinum(kk) > 0
    a1=hifreq(histart(kk));
    a2= hifreq(histop(kk));
    fprintf(1,'kk = %3i start, stop = %10.6f %10.6f \n',kk,a1,a2);
    end
  kk=3;
  if hinum(kk) > 0
    a1=hifreq(histart(kk));
    a2= hifreq(histop(kk));
    fprintf(1,'kk = %3i start, stop = %10.6f %10.6f \n',kk,a1,a2);
    end
  kk=5;
  if hinum(kk) > 0
    a1=hifreq(histart(kk));
    a2= hifreq(histop(kk));
    fprintf(1,'kk = %3i start, stop = %10.6f %10.6f \n',kk,a1,a2);
    end
  kk=2;
  if hinum(kk) > 0
    a1=hifreq(histart(kk));
    a2= hifreq(histop(kk));
    fprintf(1,'kk = %3i start, stop = %10.6f %10.6f \n',kk,a1,a2);
    end
  end

quiet=quiet0;


