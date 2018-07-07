function [hinum,histart,histop,outnum,outstart,outstop,theratio]=... 
          tempratio15um(band,temp,boxcar,hifreq,outwave,prb,freqq,INratio,...
                strenqt,w_forq,w_selfq,jr,stuff); 

%% the last line of arguments is dumb; just for the code to be compatible with
%% tempratio15umNEW.m

%this is for the 660cm-1  (15 um band)
%
%  ----------------|-----------------------------------------|-------------
%                  f1                                        f2
%   k=klor*ratio                 k=kfull                        k=klor*ratio
%
%

% INratio is what comes from Strow's paper ... it could be quite small
% for the 15 um band, so we have to go further away before we can safely
% say      k=klor*ratio

%this function tells you how the full/lorentz ratio varies with temp 
%the answers are 1d matrices of 3 elements each :  
%    kk=1 : f < f1    kk=2 : f > f2   -------- do full = lorentz*ratio 
%    kk=3 : (f > f1)&(f < f2)         -------- do full mixing 
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
% hifreq  = input wavenumber 

epep=abs(hifreq(length(hifreq))-hifreq(1))/length(hifreq)/100.0;
bbx=(boxcar-1)/2;
 
%set up dummy variables if necessary
if (nargin==5) 
  prb='Q'; freqq=[1]; 
elseif (nargin==6) 
  freqq=[1]; 
  end 
matr=getmatr(band,prb,freqq); 

if ((INratio <= 0.05) & (INratio >= 0.005))
  %things could get dicey so change bounds
  matr(1,1)=matr(1,1)-40;
  matr(2,1)=matr(2,1)-40;
  matr(3,1)=matr(3,1)+40;
  matr(4,1)=matr(4,1)+40;
  end

%matr(:,1)
%ZZ
if (band ~= 2350)
  matr(1,1)=min(freqq)-55;
  matr(2,1)=min(freqq)-40;
  matr(3,1)=max(freqq)+40;
  matr(4,1)=max(freqq)+55;
  end

if ((band == 2320) & (prb == 'P'))
  matr(1,1)=min(freqq)-55;
  matr(2,1)=min(freqq)-40;
  matr(3,1)=2500;
  matr(4,1)=2515;
  end

if ((band == 2321) & (prb == 'P'))
  matr(1,1)=min(freqq)-55;
  matr(2,1)=min(freqq)-40;
  matr(3,1)=2350;
  matr(4,1)=2365;
  end

%%%%%%% do the lower OUTPUT resolution stuff first %%%%%%%%%%%%%%%%%%%%%%
%now find the indices for 3 regions of interest
%this is for the OUTPUT resolution
iiLessOut=find(outwave < matr(2,1));          numLessOut=length(iiLessOut);
iiMoreOut=find(outwave >= matr(3,1));         numMoreOut=length(iiMoreOut);
iiNumOut =find((outwave >= matr(2,1)) & (outwave < matr(3,1))); 
numOut=length(iiNumOut);

for kk = 1:3
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
%now find the indices for 3 regions of interest
%this is for the HIGH resolution
iiLessHi=find(hifreq < matr(2,1));            numLessHi=length(iiLessHi);
iiMoreHi=find(hifreq >= matr(3,1));           numMoreHi=length(iiMoreHi);
iiNumHi =find((hifreq >= matr(2,1)) & (hifreq < matr(3,1))); 
numHi=length(iiNumHi);

for kk = 1:3
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
        %fprintf(1,'adjusting histart(3) from %4i down to 1 \n',histart(kk))
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
        %fprintf(1,'adjusting histart(3)\n',histart(kk))
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
      %fprintf(1,'adjusting histart(3)\n',histart(kk))
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
        %fprintf(1,'adjusthistop(3) from %4i upto len(hifreq) \n',histop(kk))
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
        %fprintf(1,'adjusting histop(3)\n',histop(kk))
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
      %fprintf(1,'adjusting histop(3)\n',histop(kk))
      histart(kk)=-1;
      histop(kk)=-1;
      hinum(kk)=0;       %NO POINTS in vector
      end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% to print out results
pprr=-1; 
if pprr > 0 
  [matr(2,1) matr(3,1)] 
 
  summ=0; 
  for ii=1:3      
    if (outnum(ii) > 0) 
      doh=[ii  outstart(ii) outstop(ii) ... 
           outwave(outstart(ii)) outwave(outstop(ii)) outnum(ii)]; 
      summ=summ+outnum(ii); 
      fprintf(1,'%3i %6i %6i %9.7e  %9.7e  %6i \n',doh); 
      end; 
    end 
  [summ length(outwave)] 
 
  summ=0; 
  for ii=1:3      
    if (hinum(ii) > 0) 
      doh=[ii   histart(ii) histop(ii) ... 
           hifreq(histart(ii)) hifreq(histop(ii)) hinum(ii)]; 
      summ=summ+hinum(ii); 
      fprintf(1,'%3i %6i %6i %9.7e  %9.7e  %6i \n',doh); 
      end; 
    end 
  [summ length(hifreq)] 
  end 

%we have done this
%
%  ----------------|-----------------------------------------|-------------
%                  f1                                        f2
%   k=klor*ratio                 k=kfull                        k=klor*ratio
%   I                             III                                II
%
%

%want to do this similar to tempratio15unNEW.m
%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%

hinum    = [  hinum      0  0];
histart  = [  histart   -1 -1];
histop   = [  histop    -1 -1];
outnum   = [  outnum     0  0];
outstart = [  outstart  -1 -1];
outstop  = [  outstop   -1 -1];
theratio = [  theratio  -1 -1];

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
