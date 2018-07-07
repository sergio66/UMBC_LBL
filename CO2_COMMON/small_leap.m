function [weightJJ1] = small_leap(region,spread,ffff,weightJJ0);

%%smoothly makes a weight fcn from 0 to 1, bewteen x endpts given by region
%%make the transition over "spread" cm-1

  weightJJ1 = weightJJ0;
  
  dv = spread/2;

  rA(1) = region(1)+dv;  rA(2) = region(2)-dv;
  ii = find((ffff >= rA(1)) & (ffff <= rA(2)));
  weightJJ1(ii)  = 1.0;

  x1 = region(1)-dv;  x2 = region(1)+dv;
  m = (1-0)/spread; c = 0 - m*x1;
  ii = find((ffff >= x1) & (ffff <= x2));
  weightJJ1(ii)  = m*ffff(ii) + c;

  x1 = region(2)-dv;  x2 = region(2)+dv;
  m = (0-1)/spread; c = 1 - m*x1;
  ii = find((ffff >= x1) & (ffff <= x2));
  weightJJ1(ii)  = m*ffff(ii) + c;

