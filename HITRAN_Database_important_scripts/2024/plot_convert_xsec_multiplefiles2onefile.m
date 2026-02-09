addpath /home/sergio/SPECTRA/READ_XSEC
%addpath /asl/matlib/aslutil
addpath /home/sergio/git/matlabcode/matlibSergio/matlib2025/aslutil

fmin = 605; fmax = 2830;
iStart = 51;
iStop  = 81;

iStart = input('Start plotting from gas iS (iS >= 51) : ');;
iStop  = input('Stop plotting  at   gas iE (iE <= 81) : ');;

Y1 = 2012; Y2 = 2016;
Y1 = 2020; Y2 = 2024;

str1 = ['H' num2str(Y1)];
str2 = ['H' num2str(Y2)];

for gg = iStart : iStop
  %[iYes12,gf12,xhead12] = findxsec_plot(fmin,fmax,gg,2012); 
  %[iYes16,gf16,xhead16] = findxsec_plot(fmin,fmax,gg,2016);
  [iYes12,gf12,xhead12] = findxsec_plot(fmin,fmax,gg,Y1); 
  [iYes16,gf16,xhead16] = findxsec_plot(fmin,fmax,gg,Y2);


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
  for iCnt = 1 : length(xhead12.v1)
    v = xhead12.v{iCnt};
    absc = xhead12.absc{iCnt};
    semilogy(v,absc,'b.');
    absc1 = [absc1; absc];            
    v1 = [v1 v];
    p1 = [p1 xhead12.pres(iCnt)];
    t1 = [t1 xhead12.temp(iCnt)];    
    n1 = [n1 xhead12.npts(iCnt)];
    if iCnt == 1
      ind1 = [1 length(v1)];
    else
      ind1 = [ind1 length(v1)];
    end    
    hold on
  end
  for iCnt = 1 : length(xhead16.v1)
    v = xhead16.v{iCnt};
    absc = xhead16.absc{iCnt};
    semilogy(v,absc,'r');
    absc2 = [absc2; absc];
    v2 = [v2 v];
    p2 = [p2 xhead16.pres(iCnt)];
    t2 = [t2 xhead16.temp(iCnt)];    
    n2 = [n2 xhead16.npts(iCnt)];
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
  fprintf(1,'H2012 had %6i points H2016 had %6i points out of which %6i are common \n',length(v1),length(v2),length(I))
  plot(v1(i1),absc2(i2) ./ absc1(i1));
  title([num2str(gg) ' ' str2 '/' str1]);
  axis([min(v1(i1)) max(v1(i1)) 0.85 1.15]);

  figure(3); clf
  dn = 0.5:0.01:1.5;
  plot(dn,histc(absc2(i2) ./ absc1(i1),dn));
  title([num2str(gg) ' ' str2 '/' str1 ' histogram']);
  
%  figure(2); clf
%  iCnt16 = length(xhead16.v1);
%    v = xhead16.v{iCnt16};
%    absc = xhead16.absc{iCnt16};
%    semilogy(v,absc,'r');
%    hold on
%  for iCnt = 55:56
%    v = xhead12.v{iCnt};
%    absc = xhead12.absc{iCnt};
%    semilogy(v,absc,'b');
%    hold on
%  end
%  title(num2str(gg))
%  axis([605 2830 1e-28 1e-16]); grid

  figure(4); clf;  plot(xhead12.v1,xhead12.v2,'b+',xhead16.v1,xhead16.v2,'ro')

  if gg < iStop
    disp('ret'); pause
  end
  
end
