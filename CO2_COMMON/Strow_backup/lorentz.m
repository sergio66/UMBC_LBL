function [y]=lorentz(v,v0,T,mass,brd)
%function [y]=lorentz(v,v0,T,mass,brd)
%this is the Lorentz lineshape 
%v = frequency array
%v0 = center freq
%T  = temperature
%mass = molecular mass (amu)
%brd = broadening

%c2=1.4387863;          %K/ cm-1  from Genln2 manual
%factor=v.*tanh(c2*v/2/T)/(v0*tanh(c2*v0/2/T));
%factor=tanh(c2*v/2/T)/(tanh(c2*v0/2/T));
%y=g0*y1';

%lorentz
%1/pi = 0.318...
y=0.31830988*brd./(brd*brd + (v-v0).^2);
