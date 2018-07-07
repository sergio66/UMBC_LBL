function [thetotal,number]=sortwings(line,fmin,fmax,xfar)
%function [thetotal,number]=sortwings(line,fmin,fmax,xfar)
%now sort the lines into near and far wing lines

linct=line.linct;
spectra=line.wnum(1:linct);
near=((spectra >= (fmin-xfar)) & (spectra <= (fmax+xfar)));
far= ((spectra < (fmin-xfar)) |  (spectra  > (fmax+xfar)));

num_near=sum(near); num_far=sum(far);
if ((num_near+num_far) ~= linct) 
  error('the sums no jive!')
  end

thetotal=sum(near);

doh=find((spectra >= (fmin-xfar)) & (spectra <= (fmax+xfar)));
total.linct   = sum(near);
total.wnum   = line.wnum(doh);
total.stren  = line.stren(doh);
total.tprob  = line.tprob(doh);
total.abroad = line.abroad(doh);
total.sbroad = line.sbroad(doh);
total.els    = line.els(doh);
total.abcoef = line.abcoef(doh);
total.tsp    = line.tsp(doh);
total.iusgq  = line.iusgq(doh);
total.ilsgq  = line.ilsgq(doh);
total.iso    = line.iso(doh);


%now go ahead and find the number of lines in a certain strength region
min=1e-45;
max=1e-27;
ii=1;
thestr=total.stren;
okay=((thestr >= min) & (thestr <= max));
number(ii)=sum(okay);

min=1e-27;
max=1e-26;
ii=2;
while (max*1.1 < 1e-15)
  thestr=total.stren;
  okay=((thestr >= min) & (thestr <= max));
  number(ii)=sum(okay);
  ii=ii+1;
  min=min*10;
  max=max*10;
  end

min=1e-16
max=1.0
thestr=total.stren;
okay=((thestr >= min) & (thestr <= max));
number(ii)=sum(okay);

%fprintf(1,'\n number of near, far lines = %3i %3i',sum(near),sum(far));
