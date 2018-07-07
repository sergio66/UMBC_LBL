function [hinum,histart,histop,outnum,outstart,outstop,theratio]=... 
          tempratio15um(band,temp,boxcar,hifreq,outwave,prb,freqr,INratio,...
                strenqt,w_forq,w_selfq,jr,stuff);

%%%%%%% WARNING : IF YOU MAKE CHANGES HERE, MAKE CHANGES IN BLENDER15UM.M 
 
%this is for the 700cm-1  (15 um band)
%note that it is ALSO used for the strong 667 band!!!!!!!! and the other bands!
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
%
% this is new (differs from tempratio15um)
%
%    kk=4 : f1-1 < f < f1  -------- do cousin/full mixing blend
%    kk=5 : f2 < f < f2+1  -------- do cousin/full mixing blend
%
% for wavenumbers > 2421, do cousin
% for 2420 > wavenumbers > 2421, do cousin/mixing 

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
% ratio = full/lor for following 7 temperatures temp(1):temp(7) =  ... 
%                                    150 200 .. 400 
% actually, ratio is usually set to 1.0, as it is accurately computed
% using Strows eqn

% hifreq  = input wavenumber 

global quiet 

epep=abs(hifreq(length(hifreq))-hifreq(1))/length(hifreq)/100.0;
bbx=(boxcar-1)/2;
 
%set up dummy variables if necessary
%if (nargin==5) 
%  prb='Q'; freqr=[1]; 
%elseif (nargin==6) 
%  freqr=[1]; 
%  end 

%[a,b]=quicklorentz(temp,freqr,strenqt,w_forq,w_selfq,stuff);

matr=getmatr(band,prb,freqr); 

if (prb == 'p')
  prb='P';
elseif (prb == 'r')
  prb='R';
  end

%%%this is the new code, where we turn over to Cousin after 15 wavenumbers
%%% and close to bands, we have line mixing

matr(1,1)=min(freqr)-20;
matr(2,1)=min(freqr)-15;
matr(3,1)=max(freqr)+15;
matr(4,1)=max(freqr)+20;

if ((band==667)|(band==668))
%(prb ~= 'q') & (prb ~= 'Q'))
  matr(1,1)=min(freqr)-10; 
  matr(2,1)=min(freqr)-5;
  matr(3,1)=max(freqr)+5;
  matr(4,1)=max(freqr)+10;

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now check that matr(1,1) < matr(2,1)
if matr(1,1) > matr(2,1)
  tempppp=matr(2,1);
  matr(2,1)=matr(1,1);
  matr(1,1)=tempppp;
  end

%now check that matr(3,1) < matr(4,1)
if matr(4,1) > matr(4,1)
  tempppp=matr(4,1);
  matr(4,1)=matr(3,1);
  matr(3,1)=tempppp;
  end

%%%%%%% do the lower OUTPUT resolution stuff first %%%%%%%%%%%%%%%%%%%%%%
%now find the indices for 5 regions of interest
%this is for the OUTPUT resolution
iiLessOut=find(outwave < matr(2,1));          
  numLessOut=length(iiLessOut);
iiMoreOut=find(outwave >= matr(3,1));         
  numMoreOut=length(iiMoreOut);
iiNumOut =find((outwave >= matr(2,1)+1) & (outwave < matr(3,1)-1)); 
  numOut=length(iiNumOut);
iiLessMidOut =find((outwave >= matr(2,1)) & (outwave < matr(2,1)+1)); 
  numLessMidOut=length(iiLessMidOut);
iiMoreMidOut =find((outwave >= matr(3,1)-1) & (outwave < matr(3,1))); 
  numMoreMidOut=length(iiMoreMidOut);

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
iiLessHi=find(hifreq < matr(2,1));          
  numLessHi=length(iiLessHi);
iiMoreHi=find(hifreq >= matr(3,1));         
  numMoreHi=length(iiMoreHi);
iiNumHi =find((hifreq >= matr(2,1)+1) & (hifreq < matr(3,1)-1)); 
  numHi=length(iiNumHi);
iiLessMidHi =find((hifreq >= matr(2,1)) & (hifreq < matr(2,1)+1)); 
  numLessMidHi=length(iiLessMidHi);
iiMoreMidHi =find((hifreq >= matr(3,1)-1) & (hifreq < matr(3,1))); 
  numMoreMidHi=length(iiMoreMidHi);

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
        okay = -1;
        histop(kk) = histop(kk)- boxcar;
        hinum(kk)  = hinum(kk) - boxcar;
        while okay < 0
          if (histop(kk)-(boxcar-1)) > 0
            okay = +1;
          else
            histop(kk)=histop(kk)+1;
            hinum(kk)=hinum(kk)+1;
            end
          end   
        while (truth < 0)
          %[kk histop(kk) boxcar histop(kk)-(boxcar-1)]
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
pprr=quiet;
if pprr > 0
  fprintf(1,'printing out stuff in tempratio15um \n');

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

iPrint = -1;
if iPrint > 0
  hinum
  histart
  histop
  outnum
  outstart
  outstop
  theratio
  error('oof');
  end