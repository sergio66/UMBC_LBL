function scum=ral_jjohn_blend3(band,prb,ffff,scumRAL,scumJJ1,scumJJ2,scumJJ3)

%%% this function simply blends the scum(from RAL) and scumJJX (from JJOHNS)
%%% where strenr > 2.0e-19, use scumRAL (this is in CO2_MATFILES/hit2350,2351)
%%% where strenr < 2.0e-19, use scumJJ  (this is in CO2_MATFILES/hit2350,2351)

%%% band,prb  = 2350, 'R'   for example
%%% fr,strenr = line center, temperature adjusted line strength
%%% scumRAL   = linemix computations using parameters from RAL
%%% scumJJ1    = linemix computations using parameters from JJ, region 1
%%% scumJJ2    = linemix computations using parameters from JJ, region 2

    %%% R JULY 17, 2002
      %%% region == 1 ==> 2380.0-2391.0  dv = 0.5
      %%% region == 2 ==> 2391.0-2397.5  dv = 4.0
      %%% region == 3 ==> 2397.5-2405.0  dv = 4.0

    %%% P JULY 11, 2002
      %%% region == 3 ==> 2229-2239
      %%% region == 2 ==> 2239-2289
      %%% region == 1 ==> 2289-2306

scum = scumRAL;

if ((band ~= 2350) & (band ~= 2320) & (band ~= 2310))
  error('ral_jjohn_blend3 only handles bands == 2310,2320,2350!!!!');
  end

if ((prb == 'r') | (prb == 'R'))
  if (band == 2350)
    dv2 = 0.5;  region = [2380 2391];         %%july 17, 2002, damn good!
    weightJJ  = zeros(size(ffff));
    weightJJ = small_leap(region,dv2,ffff,weightJJ);
    weightRAL = 1 - weightJJ;    
    scum = scumRAL.* weightRAL + scumJJ1.* weightJJ;

    dv2 = 3.0; region = [2391 2399];         %%july 17, 2002, damn good!
    weightJJ  = zeros(size(ffff));
    weightJJ = small_leap(region,dv2,ffff,weightJJ);
    weightRAL = 1 - weightJJ;
    scum = scum.* weightRAL + scumJJ2.* weightJJ;

    dv2 = 4.0; region = [2399 2405];         %%july 17, 2002, damn good
    weightJJ  = zeros(size(ffff));
    weightJJ = small_leap(region,dv2,ffff,weightJJ);
    weightRAL = 1 - weightJJ;
    scum = scum.* weightRAL + scumJJ3.* weightJJ;
    end
  end

%pp = [exp(-scumRAL);exp(-scumJJ1);exp(-scumJJ2);exp(-scumJJ3)];
%plot(ffff,pp,ffff,exp(-scum),'.')
%blowup
%title('do you like the blending??');
%pause