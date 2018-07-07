function [y]=vhhOPTDEP(v0,T,m,brd)
%function [y]=vhhOPTDEP(v0,T,m,brd)
%this is the Van Vleck - Huber lineshape at v=v0: using the polynomial approx
%this is a vectorised version
%v0   = center freq
%T    = temperature
%mass = molecular mass (amu)
%brd  = broadening

%%%%%%% do the doppler widths first 
k=1.380658e-23;
c_light=2.99792458e8;         %ms-1

%amu=1.6605402e-27;            %nucleon mass/kg
%mass=m*amu;                   %change to kg
%alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light); 
%r2=2*log(2)*k/amu; 
r2=11526.218;     

mm=m;
[xxv,yyv]=size(v0);
[xxm,yym]=size(mm);
if ((xxv ~= xxm) | (yyv ~= yym))
  mm=m';
  end
[xxv,yyv]=size(v0);
[xxm,yym]=size(mm);
if ((xxv ~= xxm) | (yyv ~= yym))
  error('mass of isotopes, wavevector have different lengths');
  end
alpha_doppler=v0./(c_light*sqrt(r2*T./mm)); 

%%%%%%%%% do the g0 factor 
%  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895 
%g0=sqrt(log(2)/pi)/alpha_doppler;
g0= 0.83255461 * 0.5641895 ./ alpha_doppler; 

%define arrays for the new Voigt fcn 
%X1= (v-v0)/alpha_doppler*0.83255461;            %X1={ x1  x2   ...  xn} 
%X2=(-v-v0)/alpha_doppler*0.83255461;            %X2={-x1 -x2   ... -xn} 
Y1=brd./alpha_doppler*0.83255461;  %Y1={y y  y  y ... y} no need to use ONES
X2=(-2*v0)./alpha_doppler*0.83255461;            %X2={-x1 -x2   ... -xn} 
X1=zeros(size(X2));

%X1=0.0    hence y1=y2
y=vhhOPTDEPcalc(X1,X2,Y1,g0);

%no need for tanh factors as v=v0)
