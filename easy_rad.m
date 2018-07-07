function [rad] = easy_rad(fr,k,T,Tsurf,emiss,theta)

surf = emiss.*ttorad(fr,Tsurf);
mu = cos(2*pi*theta/360);
radiance = surf';

for ii = 1 : length(T)
 ti = T(ii);
 radiance  = radiance.*exp(-k(ii,:)/mu) + ...
             ttorad(fr,ti)'.*(1-exp(-k(ii,:)/mu));
 rad(ii,:) = radiance;
 end
