function [y]=doppler(v,v0,T,m,brd)
%function [y]=doppler(v,v0,T,m,brd)
%this function computes the doppler lineshape for 
%wavevector   = v 
%line center  = v0 
%atmomic mass = m (amu) 
%temperature  = T  (k)
%braodening   = brd

k=1.380658e-23; 
c_light=2.99792458e8;         %ms-1 
amu=1.6605402e-27;            %nucleon mass/kg 
mass=m*amu;
alpha=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light);   %%%%HALF WIDTH!!!!!
y=1/alpha/sqrt(pi)*exp(-log(2)*((v-v0)/alpha).^2);

