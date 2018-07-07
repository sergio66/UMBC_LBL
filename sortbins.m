function [wing]=sortbins(line,fmin,fmax,xin,xout)

%function [wing]=sortbinss(line,fmin,fmax,xin,xout)
% if xin > 0.0 then output the lines that lie between
%     [fmin-xout fmin-xin] U [fmax+xin fmax+xout]
% else if xin <= 0.0 then output the lines that lie between
%     [fmin-xout fmax+xout]

linct=line.linct;
spectra=line.wnum(1:linct);
if (xin > 0.0) 
  left =find((spectra >= (fmin-xout)) & (spectra <  (fmin-xin)));
  right=find((spectra >  (fmax+xin)) &  (spectra <= (fmax+xout)));
  theindex=[left right];
else
  theindex=find((spectra >= (fmin-xout)) & (spectra <= (fmax+xout)));
  end

wing.linct  = length(theindex);
wing.wnum   = line.wnum(theindex);
wing.stren  = line.stren(theindex);
wing.tprob  = line.tprob(theindex);
wing.abroad = line.abroad(theindex);
wing.sbroad = line.sbroad(theindex);
wing.els    = line.els(theindex);
wing.abcoef = line.abcoef(theindex);
wing.tsp    = line.tsp(theindex);
wing.iusgq  = line.iusgq(theindex);
wing.ilsgq  = line.ilsgq(theindex);
wing.iso    = line.iso(theindex);

%fprintf(1,'\n number of near, far lines = %3i %3i',sum(near),sum(far));
