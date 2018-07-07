function [near_wing,far_wing]=sortwings(line,fmin,fmax,xfar)
%function [near_wing,far_wing]=sortwings(line,fmin,fmax,xfar)
%now sort the lines into near and far wing lines

linct=line.linct;
spectra=line.wnum(1:linct);
near=((spectra >= (fmin-xfar)) & (spectra <= (fmax+xfar)));
far= ((spectra < (fmin-xfar)) |  (spectra  > (fmax+xfar)));

num_near=sum(near); num_far=sum(far);
if ((num_near+num_far) ~= linct) 
  error('the sums no jive!')
  end

doh=find((spectra >= (fmin-xfar)) & (spectra <= (fmax+xfar)));
near_wing.linct   = sum(near);
near_wing.wnum   = line.wnum(doh);
near_wing.stren  = line.stren(doh);
near_wing.tprob  = line.tprob(doh);
near_wing.abroad = line.abroad(doh);
near_wing.sbroad = line.sbroad(doh);
near_wing.els    = line.els(doh);
near_wing.abcoef = line.abcoef(doh);
near_wing.tsp    = line.tsp(doh);
near_wing.iusgq  = line.iusgq(doh);
near_wing.ilsgq  = line.ilsgq(doh);
near_wing.iso    = line.iso(doh);

clear doh;
doh=find((spectra < (fmin-xfar)) |  (spectra  > (fmax+xfar)));
far_wing.linct   = sum(far);
far_wing.wnum   = line.wnum(doh);
far_wing.stren  = line.stren(doh);
far_wing.tprob  = line.tprob(doh);
far_wing.abroad = line.abroad(doh);
far_wing.sbroad = line.sbroad(doh);
far_wing.els    = line.els(doh);
far_wing.abcoef = line.abcoef(doh);
far_wing.tsp    = line.tsp(doh);
far_wing.iusgq  = line.iusgq(doh);
far_wing.ilsgq  = line.ilsgq(doh);
far_wing.iso    = line.iso(doh);

%fprintf(1,'\n number of near, far lines = %3i %3i',sum(near),sum(far));
