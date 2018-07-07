function   scum = ral_jjohn_blend2(band,ffff,scumRAL,scumJJ1,scumJJ2)

%%% this function simply blends the scum(from RAL) and scumJJX (from JJOHNS)
%%% where strenr > 3.8e-20, use scumRAL (this is in CO2_MATFILES/hit2350,2351)
%%% where strenr < 3.8e-20, use scumJJ  (this is in CO2_MATFILES/hit2350,2351)

%%% band,prb  = 2350, 'R'   for example
%%% fr,strenr = line center, temperature adjusted line strength
%%% scumRAL   = linemix computations using parameters from RAL
%%% scumJJ1    = linemix computations using parameters from JJ, region 1
%%% scumJJ2    = linemix computations using parameters from JJ, region 2

    %%% JULY 4, 2002
      %%% region == 1 ==> 2380-2391
      %%% region == 2 ==> 2391-2400

if (band ~= 2350)
  error('ral_jjohn_blend2 needs band == 2350');
else
  dv2 = 0.5;
  weightJJ  = zeros(size(ffff));
  region = [2380 2391];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  region = [2249 2294];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ1.* weightJJ;

  dv2 = 2.0;
  %dv2 = 0.5;
  weightJJ  = zeros(size(ffff));
  region = [2391 2400];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  region = [2229 2249];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scum.* weightRAL + scumJJ2.* weightJJ;
  end

%plot(ffff, [exp(-scumRAL);  exp(-scumJJ1); exp(-scumJJ2); exp(-scum)]);
%blowup
%title('do you like the blending??');
%pause
