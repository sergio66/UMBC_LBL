addpath /home/sergio/SPECTRA/READ_XSEC
addpath /asl/matlib/aslutil
fmin = 605; fmax = 2830;
iStart = 51;
iStop  = 81;

iStart = input('Start plotting from gas iS (iS >= 51) : ');;
iStop  = input('Stop plotting  at   gas iE (iE <= 81) : ');;

%% recall H2016/plot_convert_xsec_multiplefiles2onefile.m shows diff between H2012 and H2016
%% so let H2020/plot_convert_xsec_multiplefiles2onefile.m shows diff between H2016 and H2020
%%
%% [sergio@taki-usr2 H2020]$ more cper_xsec_from_H2016.sc
%% # ls XSEC       gives the new files Chris downloaded
%% # -rw-rw-r-- 1 chepplew pi_strow    894708 Nov 18 10:09 N2O5.zip             g62
%% # -rw-rw-r-- 1 chepplew pi_strow  97303574 Nov 18 10:07 CF4.zip              g54
%% # -rw-rw-r-- 1 chepplew pi_strow 252739555 Nov 18 10:04 SF6.zip              g81
%% # -rw-rw-r-- 1 chepplew pi_strow 120517994 Nov 18 09:59 HCFC-141b.zip        g67  TIS IS IERD, does not have *TORR*
%% # -rw-rw-r-- 1 chepplew pi_strow 122215923 Nov 18 09:51 CFC-11.zip           g51

for gg = iStart : iStop
  %[iYes12,gf12,xhead12] = findxsec_plot(fmin,fmax,gg,2012); 
  [iYes16,gf16,xhead16] = findxsec_plot(fmin,fmax,gg,2016); 
  [iYes20,gf20,xhead20] = findxsec_plot(fmin,fmax,gg,2020); 

  clear v1 absc1 v2 absc2
  v1 = [];
  v2 = [];
  absc1 = [];
  absc2 = [];
  p1 = [];
  p2 = [];
  t1 = [];
  t2 = [];
  n1 = [];
  n2 = [];
  ind1 = [];
  ind2 = [];
  
  figure(1); clf;
  for iCnt = 1 : length(xhead16.v1)
    v = xhead16.v{iCnt};
    absc = xhead16.absc{iCnt};
    semilogy(v,absc,'b.');
    absc1 = [absc1; absc];            
    v1 = [v1 v];
    p1 = [p1 xhead16.pres(iCnt)];
    t1 = [t1 xhead16.temp(iCnt)];    
    n1 = [n1 xhead16.npts(iCnt)];
    if iCnt == 1
      ind1 = [1 length(v1)];
    else
      ind1 = [ind1 length(v1)];
    end    
    hold on
  end
  for iCnt = 1 : length(xhead20.v1)
    v = xhead20.v{iCnt};
    absc = xhead20.absc{iCnt};
    semilogy(v,absc,'r');
    absc2 = [absc2; absc];
    v2 = [v2 v];
    p2 = [p2 xhead20.pres(iCnt)];
    t2 = [t2 xhead20.temp(iCnt)];    
    n2 = [n2 xhead20.npts(iCnt)];
    if iCnt == 1
      ind2 = [1 length(v2)];
    else
      ind2 = [ind2 length(v2)];
    end    
    hold on
  end
  hold off
  axis([605 2830 1e-28 1e-16]); grid
  axis([605 1505 1e-28 1e-16]); grid  
  title(num2str(gg))

%whos v1 v2 p1 p2 t1 t2 n1 n2 absc1 absc2 ind1 ind2
%% ahhhhh bloody points are not correctly sorted oh well
figure(2); clf; plot(p1,t1,'bo',p2,t2,'rx')
figure(3); clf; plot(p1,p2,'bo');
figure(4); clf; plot(t1,t2,'bo');

pause(0.01)

  figure(2); clf
  [I,i1,i2] = intersect(v1,v2);
  fprintf(1,'H2016 had %6i points H2020 had %6i points out of which %6i are common \n',length(v1),length(v2),length(I))
  plot(v1(i1),absc2(i2) ./ absc1(i1));
  title([num2str(gg) ' H2020/H2016']);
  axis([min(v1(i1)) max(v1(i1)) 0.85 1.15]);

  figure(3); clf
  dn = 0.5:0.01:1.5;
  plot(dn,histc(absc2(i2) ./ absc1(i1),dn));
  title([num2str(gg) ' H2020/H2016 histogram']);
  
%  figure(2); clf
%  iCnt20 = length(xhead20.v1);
%    v = xhead20.v{iCnt20};
%    absc = xhead20.absc{iCnt20};
%    semilogy(v,absc,'r');
%    hold on
%  for iCnt = 55:56
%    v = xhead16.v{iCnt};
%    absc = xhead16.absc{iCnt};
%    semilogy(v,absc,'b');
%    hold on
%  end
%  title(num2str(gg))
%  axis([605 2830 1e-28 1e-16]); grid

  figure(4); clf;  plot(xhead16.v1,xhead16.v2,'b+',xhead20.v1,xhead20.v2,'ro')

  if gg < iStop
    disp('ret'); pause
  end
  
end
