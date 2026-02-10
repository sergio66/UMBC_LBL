gf = 53;
cd /asl/data/hitran/HITRAN2008/IR-XSect/Compressed-files/Junk/
cd /asl/data/hitran/HITRAN08_SERGIO/Xsec/

gd = [hitranpath '/H2024/IR-XSect/Compressed-files/'];              idd = 2024;
cder = ['cd ' gd];
eval(cder)

%fname = 'CFC-13_IR01.xsc'; [d53,i53] = info_read_xsec(fname);

gf = fname;
cd /home/sergio/HITRAN2UMBCLBL/READ_XSEC
ff = 605 : 25 : 1555;

for ii = 1 : length(ff)-1
  v1 = ff(ii); v2 = ff(ii+1); dv = 1; tp = 300; pL = 1000; db = 1;
  [absc, vgrid] = calc_xsec(gf, v1, v2, dv, tp, pL, db);
  plot(vgrid,absc); 
  pause
  end
