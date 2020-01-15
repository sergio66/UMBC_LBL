iCO2 = +1;
if iCO2 > 0
  iaBand = [618 648 662 667 720 791 1932 2080 2129 ...
            668 740 2093 ...
   	    2350 2351 2352 2354 ...
	    2320 2321 2322 ...
	    2310 2311];

  exchangelinecenters = -1;

  vers = 2012;
  vers = 2016;

  strvers = num2str(vers); strvers = strvers(3:4);

  for ii = 1 : length(iaBand)
    band = iaBand(ii);
    strengthM = 0.0;
    homepath = '/home/sergio/SPECTRA/';
    hitlin_fname = ['/asl/data/hitran/h' strvers '.by.gas/g2.dat'];
    xfar = 500;
    fprintf(1,'doing band %4i \n',band)
    [line1,line1A] = co2lines(band,['h' strvers]);
    [PQRline] = makeDAVEhitlin(band,vers,strengthM,homepath,hitlin_fname,xfar,exchangelinecenters);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iCO = -1;
if iCO > 0
  iaBand = [2150];

  exchangelinecenters = -1;

  vers = 2012;
  vers = 2016;

  strvers = num2str(vers); strvers = strvers(3:4);

  for ii = 1 : length(iaBand)
    band = iaBand(ii);
    strengthM = 0.0;
    homepath = '/home/sergio/SPECTRA/';
    hitlin_fname = ['/asl/data/hitran/h' strvers '.by.gas/g5.dat'];
    xfar = 500;
    fprintf(1,'doing band %4i \n',band)
    [line1,line1A] = colines(band,['h' strvers]);
    [PQRline] = makeCOhitlin(band,vers,strengthM,homepath,hitlin_fname,xfar);
  end
end
