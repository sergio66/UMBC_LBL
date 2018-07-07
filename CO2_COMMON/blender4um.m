function [a,b]=blender4um(band,prb,freqr,jr,region);

%works directly as in tempratio4um.m, by calling find_regions

%%given band and prb, this function returns
%%  (a) the wavenumber where we stop going full mix (f1 for 'P' and fr for 'R')
%%  (b) the width of the blending region
%
%  -------|------------|----------------------------|-------------|----------
%        f1-1          f1                          f2           f2+1
%  
%  cousin   full/cousin          full               full/cousin     cousin    
%
%   I          IV                 III                   V            II
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  iiLessOut=find(outwave < matr(2,1));                               %%Reg 1 
%    numLessOut=length(iiLessOut); 
%  iiMoreOut=find(outwave >= matr(3,1));                               %%Reg 2 
%   numMoreOut=length(iiMoreOut); 
%  iiNumOut =find((outwave >= matr(2,1)+1) & (outwave < matr(3,1)-1)); %%Reg 3 
%  iiNumOut =find((outwave >= matr(2,1)+20) &(outwave < matr(3,1)-20));%%Reg 3 
%    numOut=length(iiNumOut); 
%  iiLessMidOut =find((outwave >= matr(2,1)) &(outwave < matr(2,1)+1));%%Reg 4 
%    numLessMidOut=length(iiLessMidOut);  
%  iiMoreMidOut =find((outwave >= matr(3,1)-1) &(outwave < matr(3,1)));%%Reg 5 
%    numMoreMidOut=length(iiMoreMidOut); 

find_regions;

%%set up the width of the blending regions
if band ~= 2350
  b = 1.0;
else
  b = 20.0;
  end

%%set up the wavenumber where we stop doing fullmixing
if (region == 4)
  a=matr(2,1)+b;
elseif (region == 5)
  a=matr(3,1)-b;
else
  error('in blender 4um.m .... need region = 4 or 5 \n');
  end

