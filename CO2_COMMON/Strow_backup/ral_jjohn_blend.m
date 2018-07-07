function   scum = ral_jjohn_blend(band,prb,fr,strenrt,ffff,scumRAL,scumJJ)

%%% this function simply blends the scum(from RAL) and scumJJ (from JJOHNS)
%%% where strenr > 3.8e-20, use scumRAL (this is in CO2_MATFILES/hit2350,2351)
%%% where strenr < 3.8e-20, use scumJJ  (this is in CO2_MATFILES/hit2350,2351)

%%% band,prb  = 2350, 'R'   for example
%%% fr,strenr = line center, temperature adjusted line strength
%%% ffff      = kcarta output frequency (before boxcar integration)
%%% scumRAL   = linemix computations using parameters from RAL
%%% scumJJ    = linemix computations using parameters from JJ

scum = scumRAL;

if ((band ~= 2351) & (band ~= 2321) & (band ~= 2311))
  error('ral_jjohn_blend only handles bands == 2311,2321,2351!!!!');

elseif ((band == 2351) & ((prb == 'p') | (prb == 'P')))
  weightJJ  = zeros(size(ffff));
  dv2 = 4.0; region = [2270 2295];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ.* weightJJ;
elseif ((band == 2351) & ((prb == 'r') | (prb == 'R')))
  fprintf(1,'band prb = %5i %s \n',band,prb);
  weightJJ  = zeros(size(ffff));
  dv2 = 4.0; region = [2255 2295];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ.* weightJJ;

elseif ((band == 2321) & ((prb == 'p') | (prb == 'P')))
  weightJJ  = zeros(size(ffff));
  dv2 = 4.0; region = [2270 2275];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ.* weightJJ;
elseif ((band == 2321) & ((prb == 'r') | (prb == 'R')))
  scum = scumRAL;              %%%very gradual change in intensity

elseif ((band == 2311) & ((prb == 'p') | (prb == 'P')))
  weightJJ  = zeros(size(ffff));
  dv2 = 4.0; region = [2258 2265];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ.* weightJJ;
elseif ((band == 2311) & ((prb == 'r') | (prb == 'R')))
  weightJJ  = zeros(size(ffff));
  dv2 = 4.0; region = [2255 2265];
  weightJJ = small_leap(region,dv2,ffff,weightJJ);
  weightRAL = 1 - weightJJ;
  scum = scumRAL.* weightRAL + scumJJ.* weightJJ;

  end

%plot(ffff, exp(-scumRAL),'.', ffff, exp(-scumJJ),'.', ffff, exp(-scum));
%blowup
%title('do you like the blending??');
%pause