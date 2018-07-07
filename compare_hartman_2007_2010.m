dd = load('IPFILES/co2two');

load hartnew_compare_2007_2010

iDo = -1;
if iDo > 0
  clear topts
  topts.ffin = 0.0025;
  topts.nbox = 1;
  topts.iVersHartmann = 2007;
  [w,d07A] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoA',topts);
  [w,d07B] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoB',topts);
  save hartnew_compare_2007_2010 w d07*
  end

iDo = -1;
if iDo > 0
  clear topts
  topts.ffin = 0.0025;
  topts.nbox = 1;
  topts.iVersHartmann = 2010;
  topts.WV_partialpressure = 0.0*dd(1,2);
  topts
  [w,d10A] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoA',topts);
  topts.WV_partialpressure = 0.0*dd(2,2);
  [w,d10B] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoB',topts);
  save hartnew_compare_2007_2010 w d07* d10*
  end

iDo = -1;
if iDo > 0
  clear topts
  topts.ffin = 0.0025;
  topts.nbox = 1;
  topts.iVersHartmann = 2010;
  topts.WV_partialpressure = 0.2*dd(1,2);
  topts
  [w,dx10A] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoA',topts);
  topts.WV_partialpressure = 0.2*dd(2,2);
  [w,dx10B] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoB',topts);
  save hartnew_compare_2007_2010 w d07* d10* dx10*
  end

iDo = +1;
if iDo > 0
  clear topts
  topts.mainloop = -1;     %% this way, ignore my "remaining lines" loop
  topts.ffin = 0.0025;
  topts.nbox = 1;
  topts.iVersHartmann = 2010;
  topts.WV_partialpressure = 0.2*dd(1,2);
  topts
  [w,dy10A] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoA',topts);
  topts.WV_partialpressure = 0.2*dd(2,2);
  [w,dy10B] =run8co2_hartmann(2,605,2830,'IPFILES/co2twoB',topts);
  save hartnew_compare_2007_2010 w d07* d10* dx10* dy10*
  end

plot(w,d10A./d07A,w,dx10A./d07A)
plot(w,[d07A; d10A; dx10A]); grid