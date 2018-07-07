function [a,b]=blender15um(band,prb,freqr,region);
%%copied directly from tempratio15um.m

%%%%%%% WARNING : IF YOU MAKE CHANGES HERE, MAKE CHANGES IN TEMPRATIO15UM.M

%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  iiLessOut=find(outwave < matr(2,1));                                  %%Reg 1 
%    numLessOut=length(iiLessOut); 
%  iiMoreOut=find(outwave >= matr(3,1));                                 %%Reg 2 
%   numMoreOut=length(iiMoreOut); 
%  iiNumOut =find((outwave >= matr(2,1)+1) & (outwave < matr(3,1)-1));     %%Reg 3 
%    numOut=length(iiNumOut); 
%  iiLessMidOut =find((outwave >= matr(2,1)) & (outwave < matr(2,1)+1)); %%Reg 4 
%    numLessMidOut=length(iiLessMidOut);  
%  iiMoreMidOut =find((outwave >= matr(3,1)-1) & (outwave < matr(3,1))); %%Reg 5 
%    numMoreMidOut=length(iiMoreMidOut); 

%%set up the width of the blending regions
b=1.0;

%%set up the wavenumber where we stop doing fullmixing
if (region == 4)
  a=matr(2,1)+b;
elseif (region == 5)
  a=matr(3,1)-b;
else
  error('in blender 15um.m .... need region = 4 or 5 \n');
  end

