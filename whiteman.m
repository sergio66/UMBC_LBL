f0 = [351 354.7 371 375.4 382 387 403 407.5  532.4 580.4 607.8 660.9]; %%%in nm
f1 = f0 + 0.3;                                    %%%wavelength interval

f0 = f0 * 1e-9 * 1e6;                             %%% change to um
f1 = f1 * 1e-9 * 1e6;                             %%% change to um

f0 = 10000./f0;   %%%change to cm-1
f1 = 10000./f1;   %%%change to cm-1

f0 = ceil(f0);
f1 = floor(f1);

for ii = 1 : length(f0)
  figure(1)
  fa = f1(ii);
  fb = f0(ii);

  profname = 'IPFILES/std_water';
  clear topts;
  topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
  hitranpath.m
  topts.HITRAN = HITRAN;  
  [w,k1line] = run7water(1,fa,fb,profname,topts);
  [w,k1con]  = run7watercontinuum(1,fa,fb,profname);
  k1 = k1line + k1con;

  profname = 'IPFILES/std_co2';
  [w,k2] = run7(2,fa,fb,profname);

  profname = 'IPFILES/std_ozone';
  [w,k3] = run7(3,fa,fb,profname);

  profname = 'IPFILES/std_n2o';
  [w,k4] = run7(4,fa,fb,profname);

  profname = 'IPFILES/std_co';
  [w,k5] = run7(5,fa,fb,profname);

  profname = 'IPFILES/std_ch4';
  [w,k6] = run7(6,fa,fb,profname);

  clear topts; topts.CKD = +1;
  profname = 'IPFILES/std_o2';
  [w,k7] = run7(7,fa,fb,profname,topts);

  clear topts; topts.CKD = +1;
  profname = 'IPFILES/std_n2';
  [w,k22] = run7(22,fa,fb,profname,topts);

  ktotal = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k22;

  figure(2)
  opticaldepth = sum(ktotal);
  plot(10000./w * 1000,exp(-opticaldepth)); pause(1)

  end
