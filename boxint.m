function [z]=boxint(y,nbox);
%function [z]=boxint(y,nbox); does a boxcar integration

z=zeros(1,length(y)/nbox);
zlen=length(z);
for jj=1:zlen
  ind=(1:nbox)+(jj-1)*nbox;
  z(jj)=sum(y(ind));
  end
z=z/nbox;

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
