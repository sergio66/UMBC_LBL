function [y]=vhh1(v,v0,T,m,brd)
%function [y]=vhh1(v,v0,T,m,brd)
%this is the Van Vleck - Huber lineshape  : uses Dave Tobin's voigt profile
%v   = frequency array
%v0  = center freq
%T   = temperature
%m   = molecular mass (amu)
%brd = lorentz broadening

%%%%%%% do the doppler widths first
k=1.380658e-23; 
c_light=2.99792458e8;         %ms-1 
amu=1.6605402e-27;            %nucleon mass/kg 
mass=m*amu;

%alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light);
%r2=2*log(2)*k/amu;
r2=11526.218;             
alpha_doppler=v0/c_light*sqrt(r2*T/m);

%%%%%%%%% do the g0 factor
%  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
%g0=sqrt(log(2)/pi)/alpha_doppler;
g0= 0.83255461 * 0.5641895 / alpha_doppler;

%define arrays for the new Voigt fcn
X1= (v-v0)/alpha_doppler*0.83255461;            %X1={x1 x2 x3  ... xn}
X2=(-v-v0)/alpha_doppler*0.83255461;            %X1={x1 x2 x3  ... xn}
Y1=ones(size(X1))*brd/alpha_doppler*0.83255461;
                                                  %Y1={y y  y  y ... y}
y1=real(voigt1calc(X1,Y1,g0));
y2=real(voigt1calc(X2,Y1,g0));

c2=1.4387863;          %K/ cm-1  from Genln2 manual

factor=v.*tanh(c2*v/2/T)/(v0*tanh(c2*v0/2/T));
y=factor.*(y1+y2)';

