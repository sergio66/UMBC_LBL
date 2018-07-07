function [z] = boxint_many(y,nbox);

%function [z] = boxint(y,nbox); does a boxcar integration of a matrix
% nbox is number of points over which to do the boxcar

[m,n] = size(y);
if m > n
  y = y';
  end

[m,n] = size(y);

% fprintf(1,'doing boxcar over %4i points number of wavenumber points = %5i \n',nbox,n)

z     = zeros(m,n/nbox);
zlen  = n/nbox;
for jj = 1:zlen
  ind    = (1:nbox)+(jj-1)*nbox;
  summer = y(:,ind);   
  summer = sum(summer');
  z(:,jj) = summer';
  end
z = z/nbox;

%this is if we want to do trapezoid integration
%z=zeros(1,length(y)/nbox);
%zlen=length(z);
%for jj=1:zlen
%  ind=(1:nbox)+(jj-1)*nbox;
%  indend=[ind(1) ind(nbox)];
%  indmid=ind(2:nbox-1);
%  z(jj)=sum(0.5*y(indend))+sum(y(indmid));
%  end
%z=z/4.0;
