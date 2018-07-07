function [voi]=smooth_full_ratio(outwave,voi,df,outnum,outstart,outstop); 
%this functions smooths the transition regimes from k= kfull to k=klor*ratio
%it just does a linear interpolation from
%y(-df) to y(+df) where -df,...,df are the points on either side of transition
%such that they span about 0.1 cm-1

global quiet

voiorig=voi;

%outnum
%outstart
%outstop

if ((outnum(4) ~= 0)  & (outnum(3) ~= 0))
  ii=outstop(4);      %this is where we have to smooth
  if ii > df
    ind=ii-df:ii+df;
    slope=(voi(ind(length(ind)))-voi(ind(1)))/(2*df+1);
    jj=-df:df;
    voi(ind)=slope*(jj+df)+voi(ind(1));
    end
  end
if ((outnum(5) ~= 0)  & (outnum(3) ~= 0))
  ii=outstart(5);     %this is where we have to smooth
  if (ii < (length(voi)-df))
    ind=ii-df:ii+df;
    slope=(voi(ind(length(ind)))-voi(ind(1)))/(2*df+1);
    jj=-df:df;
    voi(ind)=slope*(jj+df)+voi(ind(1));
    end
  end


if ((outnum(1) ~= 0)  & (outnum(4) ~= 0))
  ii=outstop(1);      %this is where we have to smooth
  if ii > df
    ind=ii-df:ii+df;
    slope=(voi(ind(length(ind)))-voi(ind(1)))/(2*df+1);
    jj=-df:df;
    voi(ind)=slope*(jj+df)+voi(ind(1));
    end
  end
if ((outnum(5) ~= 0)  & (outnum(2) ~= 0))
  ii=outstart(2);     %this is where we have to smooth
  if ((ii < (length(voi)-df)) & (ii > df))
    ind=ii-df:ii+df;
    slope=(voi(ind(length(ind)))-voi(ind(1)))/(2*df+1);
    jj=-df:df;
    voi(ind)=slope*(jj+df)+voi(ind(1));
    end
  end

if quiet > 0
  plot(outwave,exp(-voi),outwave,exp(-voiorig))
  title('after smoothing')
  pause(0.05)
  end
