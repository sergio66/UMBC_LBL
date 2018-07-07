  iaBand = [618 648 662 667 720 791 1932 2080 2129 ...
            668 740 2093 ...
            2350 2351 2352 2354 ...
            2320 2321 2322 ...
            2310 2311];


for ii = 1 : length(iaBand)
  band = iaBand(ii);
  fname0 = ['/home/sergio/SPECTRA/CO2_MATFILES/H16/hit' num2str(band) '.mat'];
  fname1 = ['/home/sergio/SPECTRA/xtest_CO2_MATFILES/H16/hit' num2str(band) '.mat'];
  if exist(fname0) & exist(fname1)
    wah0 = load(fname0);
    wah1 = load(fname1);
    fprintf(1,'band %4i len(wah0) %5i len(wah0) %5i \n',band,length(wah0.freq),length(wah1.freq))
    fprintf(1,'  diff(stren) diff(wnum) %8.6e %8.6e \n',sum(wah0.stren-wah1.stren),sum(wah0.freq-wah1.freq))
    semilogy(wah0.freq,wah0.stren,'b.',wah1.freq,wah1.stren,'ro')
    title(num2str(band))
    disp(' ')
    pause(1)
  end
end
