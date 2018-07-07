file1a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/my049202.txt';
file2a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl069208.txt';
file3a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl079201.txt';
file4a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl079203.txt';
file5a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl089201.txt';
file6a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl109203.txt';
file7a = '/carrot/s2/strow/Tobin/Tobin_home/Co2pr/Johns_jul92/jl139201.txt';

d1a = load(file1a);
d2a = load(file2a);
d3a = load(file3a);
d4a = load(file4a);
d5a = load(file5a);
d6a = load(file6a);
d7a = load(file7a);

%%%now we wanna check how the calcs measure upto the obs
load /asl/data/ral/Co2a/file1n2_johnjohn_all.mat
load /asl/data/ral/Co2a/johnjohn_c02_k4.mat
ksim1 = kn2 + k4;                 %%set amount of CO2
ksim2 = kn2 + 1.1*k4;             %%perturbed by 1.1

%%%at ii = 5, do  [w0,y0]=ginput; to get the wavenumbers in between lines
%%%save /home/sergio/SPECTRA/CO2_RAL_FITS/inbetween4um.mat w0
load /home/sergio/SPECTRA/CO2_RAL_FITS/inbetween4um.mat
load /home/sergio/AIRSPRODUCTS/FIRSTLIGHT/chi_co2.mat
ipfile = load ('/home/sergio/SPECTRA/CO2_RAL_FITS/IPFILES/davethesis_co2');
p = ipfile(:,2)*1013; ps = ipfile(:,3)*1013; T = ipfile(:,4);

%try a triangle chi function : 1.09 at 2380, 1.13 at 2390 and 1.05 at 2405
m1 = (1.13-1.09)/(2390-2380); c1 = 1.13-m1*2390;
x1 = 2355:0.0025:2390-0.0025; y1 = m1*x1 + c1;
m2 = (1.05-1.13)/(2405-2390); c2 = 1.05-m2*2405;
x2 = 2390:0.0025:2430-0.0025; y2 = m2*x2 + c2;
xx = [x1 x2]; yy = [y1 y2]; yy(find(yy < 1)) = 1.0;
yy0 = interp1(xx,yy,w0);

for ii = 4  : 4
  if (ii == 1)
    fdata = d1a(:,1); data = d1a(:,2);
  elseif (ii == 2)
    fdata = d2a(:,1); data = d2a(:,2);
  elseif (ii == 3)
    fdata = d3a(:,1); data = d3a(:,2);
  elseif (ii == 4)
    fdata = d4a(:,1); data = d4a(:,2);
  elseif (ii == 5)
    fdata = d5a(:,1); data = d5a(:,2);
  elseif (ii == 6)
    fdata = d6a(:,1); data = d6a(:,2);
  elseif (ii == 7)
    fdata = d7a(:,1); data = d7a(:,2);
    end

  subplot(311);
  plot(fr,exp(-ksim1(ii,:)),fr,exp(-ksim2(ii,:)),fdata,data,'r'); grid 
  axis([2380 2400 0 1]); 

  orig = interp1(fr,exp(-ksim1(ii,:)),fdata);
  new = interp1(fr,exp(-ksim2(ii,:)),fdata);
  subplot(312);
  plot(fdata,data-orig,fdata,data-new,'r'); grid
  axis([2380 2400 -0.1 0.1]); 
 
  subplot(313);
  datak = -log(data); calcsk = -log(orig);
  fratio = 2380:0.1:2400;
  ratio = datak./calcsk; 
  ratio1 = interp1(fdata,ratio,fratio); 
  ratiow0 = interp1(fdata,ratio,w0); 
  chi(ii,:) = ratio1;
  str = [num2str(ii) ' : p,ps,T = ' num2str(p(ii)) ' ' num2str(ps(ii)) ' ' ...
         num2str(T(ii))];
  fprintf(1,'%s \n',str);
  %%plot(fratio,chi(ii,:),fZ1,Z,w0,ratiow0,'r.-'); title(str); grid
  plot(fratio,chi(ii,:),w0,yy0,w0,ratiow0,'r.-'); title(str); grid
  axis([2380 2400 -0.5 2.0]); 

  kcartasim = exp(-ksim1(ii,:));
  kcartasim = interp1(fr,kcartasim,fdata);
  save strow.mat fdata kcartasim data str
  end

clear
clf
load strow.mat
plot(fdata,kcartasim,fdata,data,'r');
title(str)
plot(fdata,data-kcartasim,'r');