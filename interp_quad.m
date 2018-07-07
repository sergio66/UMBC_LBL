function yspline=interp_quad(x,y,x1)

if (length(x) ~= 3)
  error('need a quad interp ==> need 3 points');
  end

%ma=1;
%if (ma > 0)        %use matrix methods directly
  matr=[x.^2;x;ones(1,3)]';
  a=inv(matr)*y';
  yspline=a(1)*x1.^2+a(2)*x1+a(3);
%else
% t=cputime;
%  del2=x(2)-x(1);   del2del2=x(2)*x(2)-x(1)*x(1);
%  del3=x(3)-x(1);   del3del3=x(3)*x(3)-x(1)*x(1);
%  gam2=y(2)-y(1);
%  gam3=y(3)-y(1);
%  a(1)=(gam2*del3-gam3*del2)/(del2del2*del3-del3del3*del2);
%  a(2)=(gam2-a(1)*del2del2)/del2;
%  a(3)=y(1)-a(1)*x(1)*x(1)-a(2)*x(1);
%  yspline=a(1)*x1.^2+a(2)*x1+a(3);
% cputime-t 
%  end



