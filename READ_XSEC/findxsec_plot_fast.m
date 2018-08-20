function iYes = findxsec_plot_fast(wn1,wn2,bands)

%% assumes you have previously done eg
%%    bands = list_bands(gg);
%% input  
%%   wn1,wn2  = start/stop chunk you want to check
%%   numbands = number of bands
%%   bands.v1, bands.v2 = arrays containing the start/stop of bands
%%
%% output
%%   iYes = +/-1 depending on whether the wavnumber chunk lies inside a band

iYes = -1;

numbands = length(bands.v1);
if numbands == 0
  wn1
  wn2
  iYes = -1;
  disp('oh oh nothing found')
end

iCnt = 1;
while iCnt <= numbands & iYes < 0
  if wn1 <= bands.v2(iCnt) & wn2 >= bands.v1(iCnt)
    iYes = 1;
  else
    iCnt = iCnt + 1;
  end
end
