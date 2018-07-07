function y = vhh_fast(v,v0,T,mass,brd)

%%% this function does a fast vhh, using lorentz stuff in the far wings
%%% idea is from J. Lenoble's book, that in far wing, voigt ~ lorentz

%v = frequency array
%v0 = center freq
%T  = temperature
%mass = molecular mass (amu)
%brd = broadening

%c2=1.4387863;          %K/ cm-1  from Genln2 manual
%factor=v.*tanh(c2*v/2/T)/(v0*tanh(c2*v0/2/T));
%factor=tanh(c2*v/2/T)/(tanh(c2*v0/2/T));
%y=g0*y1';

%lorentz 1/pi = 0.318...
%%%y=0.31830988*brd./(brd*brd + (v-v0).^2);    assume |v+/-vo| >> brd

% alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light) 
% r2=2*log(2)*k/amu 
k=1.380658e-23;
c_light=2.99792458e8;         
amu=1.6605402e-27;            %%nucleon mass/kg 

r2=11526.218;
alpha_doppler=v0/c_light*sqrt(r2*T/mass);
repwid=0.8325546/alpha_doppler;

g0= 0.83255461 * 0.5641895 / alpha_doppler;

xp = (v-v0); yp = 0.31830988 * brd./(xp.^2 + 1e-38);
xm = (v+v0); ym = 0.31830988 * brd./(xm.^2);

%do the vhh part .. see FORTRANLINUX/vhh2RI.f

c2   = 1.4387863;        %%K/ cm-1  from Genln2 manual  
c2_0 = 0.5*c2*v0/T; 
c2_0 = v0*tanh(c2_0);  
fact = repwid*0.5641895/c2_0;
  
c2=0.5*c2/T;  
factor = fact*v.*tanh(c2*v);
y      = factor.*(yp+ym)/g0;
