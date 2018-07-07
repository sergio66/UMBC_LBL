function band = reverse_engineer_co2band(xx);

band = NaN;

isotope = xx(1);
v_u     = xx(2);
v_l     = xx(3);

if     (v_l == 2 & v_u == 3 & isotope == 1)
  band = 618;
elseif (v_l == 1 & v_u == 2 & isotope == 2)
  band = 648;
elseif (v_l == 1 & v_u == 2 & isotope == 3)
  band = 662;
elseif (v_l == 1 & v_u == 2 & isotope == 1)
  band = 667;
elseif (v_l == 2 & v_u == 5 & isotope == 1)
  band = 720;

elseif (v_l == 3 & v_u == 8 & isotope == 1)
  band = 791;
elseif (v_l == 1 & v_u == 8 & isotope == 1)
  band = 2080;

%find all PQR lines for isotope 1 : these are PQR_deltpi
elseif (v_l == 2 & v_u == 4 & isotope == 1)
   band = 668;
elseif (v_l == 4 & v_u == 8 & isotope == 1)
  band = 740;
elseif (v_l == 2 & v_u == 14 & isotope == 1)
  band = 2093;

%find all PQR lines for isotope 1 : these are PQR_sigsig
elseif (v_l == 1 & v_u == 9 & isotope == 1)
 band = 2350;
elseif (v_l == 1 & v_u == 9 & isotope == 2)
  band = 2351;
elseif (v_l == 1 & v_u == 9 & isotope == 3)
  band = 2352;
elseif (v_l == 3 & v_u == 23 & isotope == 1)
 band = 2353;
elseif (v_l == 5 & v_u == 25 & isotope == 1)
 band = 2354;

%find all PQR lines for isotope 1 : these are PQR_pipi
elseif ( v_l == 2 & v_u == 16 & isotope == 1)
  band = 2320; 
elseif ( v_l == 2 & v_u == 16 & isotope == 2)
  band = 2321; 
elseif ( v_l == 2 & v_u == 16 & isotope == 3)
  band = 2322; 

%find all PQR lines for isotope 1 : these are PQR_deltdelt
elseif ( v_l == 4 & v_u == 24 & isotope == 1)
  band = 2310; 
elseif ( v_l == 4 & v_u == 24 & isotope == 2)
  band = 2311; 
  end

%[v1 v2 L v3 r] <----- [v1 v2 L v3 r]
%class 5    CO2 
%c                                       .......1   .......2   .......3    
%      DATA (AVIB(I),I=171,229)/        '   00001','   01101','   10002', 
%c 
%c      .......4   .......5   .......6   .......7   .......8   .......9 
%     +'   02201','   10001','   11102','   03301','   11101','   00011', 
%c 
%c      ......10   ......11   ......12   ......13   ......14   ......15 
%%     +'   20003','   12202','   20002','   04401','   12201','   20001', 
%c 
%c      ......16   ......17   ......18   ......19   ......20   ......21 
%     +'   01111','   21103','   13302','   21102','   05501','   13301', 
%c 
%c      ......22   ......23   ......24   ......25   ......26   ......27 
%     +'   21101','   10012','   02211','   10011','   30004','   22203', 
%c 
%c      ......28   ......29   ......30   ......31   ......32   ......33 
%     +'   14402','   30003','   22202','   06601','   30002','   14401', 
%c 
%c      ......34   ......35   ......36   ......37   ......38   ......39 
%     +'   22201','   30001','   11112','   03311','   11111','   00021', 
%c 
%c      ......40   ......41   ......42   ......43   ......44   ......45 
%     +'   31104','   31103','   31102','   20013','   12212','   23301', 
%c 
%c      ......46   ......47   ......48   ......49   ......50   ......51 
%     +'   31101','   04411','   20012','   12211','   20011','   01121', 
%c 
%c      ......52   ......53   ......54   ......55   ......56   ......57 
%     +'   40004','   32203','   21113','   40002','   13312','   05511', 
%c 
%c      ......58   ......59 
%     +'   21112','   13311'/ 

