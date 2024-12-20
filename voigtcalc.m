function [y]=voigt1calcG(x1,y1,g0) 
%function [y]=voigt1calcF(x1,y1,g0) 
%this function computes REAL part of voigt lineshape for wavevector x, 
% real numbere y 
%using the polynomial approx to the voigt fcns 
% 
%faster version of voigtcalc.m
 
%more accurate ones found from find_root.m 
%a1= -1.35094358543273 + 0.37861161238627i;  
%a2= -1.21498213255730 + 1.23588765343593i;  
%c1 =  0.59059188440886 - 1.18584332504043i 
%c2 = -0.30849709263498 + 0.02098588080036i 
 
a=[-1.21498213255730 -1.35094358543273 -1.21498213255730 -1.35094358543273]; 
b=[ 1.23588765343593  0.37861161238627 -1.23588765343593 -0.37861161238627]; 
c=[-0.30849709263498  0.59059188440886 -0.30849709263498  0.59059188440886]; 
d=[ 0.02098588080036 -1.18584332504043 -0.02098588080036  1.18584332504043]; 
 
z=zeros(size(x1)); 
for ii=1:4 
  z=z+(c(ii)*(y1-a(ii))+d(ii)*(x1-b(ii)))./((y1-a(ii)).^2 + (x1-b(ii)).^2); 
  end 
y=g0*z; 
 
