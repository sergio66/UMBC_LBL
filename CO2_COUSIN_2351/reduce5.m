function [f2,trd2,back2,nit2,zz2]=reduce(f,trd,back,nit,zz,lowerf);

f2    = lower_tobin(f,lowerf);
trd2  = lower_tobin(trd,lowerf);
back2 = lower_tobin(back,lowerf);
nit2  = lower_tobin(nit,lowerf);
zz2   = lower_tobin(zz,lowerf);

