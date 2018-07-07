function [wr,wi]=voigt1calc(x,y,g0);
%function [wr,wi]=voigt1calc(x,y,g0);
%computes real part and imag part of voigt fcn (imag part for line mixing)
% 
%x,y have come in from the calling function
%x=(v-v0)/alpha_doppler*0.83255461;           
%y=ones(size(X1))*brd/alpha_doppler*0.83255461;  


T=[0.314240376  0.947788391   1.59768264    2.27950708     3.02063703 ...
   3.8897249];
C=[1.01172805  -0.75197147    1.2557727E-2  1.00220082E-2 -2.42068135E-4  ...
   5.00848061E-7];
S=[1.393237     0.231152406  -0.155351466   6.21836624E-3   9.19082986E-5  ...
   -6.27525958E-7];

region1= (y>0.85) | (abs(x)<(18.1*y+1.65));
region1_i=find(region1);
region2_i=find(~region1);
xr1=x(region1_i);
yr1=y(region1_i);
xr2=x(region2_i);
yr2=y(region2_i);

% Do region 1 first

if length(xr1) > 0 
   n1=length(xr1);
   wr1=zeros(1,n1);
   wi1=zeros(1,n1);
   y1=yr1+1.5;
   y2=y1.^2;
   for i=1:6 
      r=xr1-T(i);
      d=ones(1,n1)./(r.^2+y2);
      d1=y1.*d;
      d2=r.*d;
      r=xr1+T(i);
      d=ones(1,n1)./(r.^2+y2);
      d3=y1.*d;
      d4=r.*d;
      wr1=wr1+C(i)*(d1+d3)-S(i)*(d2-d4);
      wi1=wi1+C(i)*(d2+d4)+S(i)*(d1-d3);
   end
   clear r y1 y2 d d1 d2 d3 d4 
end

% Now do region 2

if length(xr2) > 0
   n2=length(xr2);
   wr2=zeros(1,n2);
   wi2=zeros(1,n2);
   
   y1=yr2+1.5;
   y2=y1.^2;
   if abs(xr2)<12 
%      plot(abs(xr2))
%      pause
      wr2=exp(-xr2.^2); end
   y3=yr2+3;
   for i=1:6 
      r=xr2-T(i);
      r2=r.^2;
      d=ones(1,n2)./(r2+y2);
      d1=y1.*d;
      d2=r.*d;
      wr2=wr2+yr2.*(C(i)*(r.*d2-1.5*d1)+S(i)*(y3.*d2))./(r2+2.25);
      r=xr2+T(i);
      r2=r.^2;
      d=ones(1,n2)./(r2+y2);
      d3=y1.*d;
      d4=r.*d;
      wr2=wr2+yr2.*(C(i)*(r.*d4-1.5*d3)-S(i)*(y3.*d4))./(r2+2.25);
      wi2=wi2+C(i)*(d2+d4)+S(i)*(d1-d3);
   end 
end
 
n=length(x);
wr=zeros(n,1);
wi=zeros(n,1);

if length(xr1) > 0
   wr(region1_i)=wr1;
   wi(region1_i)=wi1;
end
if length(xr2) > 0 
   wi(region2_i)=wi2;
   wr(region2_i)=wr2;
end

%[max(abs(wr))]
%plot(x,wr,x,wi)

wr=wr*g0;
wi=wi*g0;
