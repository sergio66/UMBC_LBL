function z = fasthankel_1(x)

%% birnbaum code uses chi1=zz.*(-pi/2).*real(besselh(1,1,sqrt(-1)*zz)).*ex;
%% so this does J1,Y1 for |x| < 8
%% see Numerical Recipes

ii = (abs(x) > 8);
if (sum(ii) > 0)
  z = real(besselh(1,1,sqrt(-1)*x));
else
  %%% j1
  y = x.*x;
  ans1 = x.*(7236261432 + y.*(-7895059235 + y.*(242396853.1 ...
        +y.*(-2972611.439+y.*(15704.48260+y*(-30.16036606))))));    
  ans2 =     144725228442.0+y.*(2300535178.0+y.*(18583304.74 ...
        +y.*(99447.43394+y.*(3769991397+y*1.0))));
  j1 = ans1./ans2;

  %%% y1
  y = x.*x;
  ans1 =  x.*(-0.4900604943e13 + y.*(0.1275274390e13 ...
         +y.*(-.05153438139e11 + y.*(0.7349264551e9  ...
         +y.*(-0.4237922726e7 +  y*0.8511937935e4)))));
  ans2 =      0.2499580570e14+y.*(0.424419664e12 ...
         +y.*(
  z = j1 + i * y1;
  z = z.*exp(-i * x);
  end
