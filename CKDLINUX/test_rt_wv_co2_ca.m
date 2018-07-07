figure(1);
f0 = 2300:2600;
f0 = 2000:2800;

CA300 = wv_co2_ca(300,1.0,f0);
iWrite = input('write out CA_WV_CO2 (-1/+1) and quit? ');
if iWrite > 0
  fdir = ['/home/sergio/SPECTRA/CKDLINUX/'];
  fname = 'ca_wv_co2_forkcarta_2000_2900.dat';
  dtype = 'ieee-le';  
  fid=fopen([fdir '/' fname], 'w' ,dtype);
  
  filemark = 4;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,length(f0),'integer*4');
  fwrite(fid,filemark,'integer*4');

  filemark = 4*length(f0);
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,f0,'real*4');
  fwrite(fid,filemark,'integer*4');

  filemark = 4*length(f0);
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,CA300,'real*4');
  fwrite(fid,filemark,'integer*4');

  fclose(fid);
  error('ok done')
end


  addpath /home/sergio/SPECTRA
  hitname = '/asl/data/hitran/h16.by.gas/g2.dat';
  [lineORIG,hitran_version,hlist_qtips] = hitread(2000,2600,0,2,hitname,-1);
  iso = find(lineORIG.iso == 1 | lineORIG.iso == 2);
  yyaxis left
  semilogy(lineORIG.wnum,lineORIG.stren,'b.',lineORIG.wnum(iso),lineORIG.stren(iso),'cx')
  yyaxis right
  semilogy(f0,CA300,'r')

figure(2)
hatran = load('co2_wv_chi_ha_tran.dat');
errorbar(hatran(:,1),hatran(:,2),hatran(:,3),'color','b','linewidth',2);
hold on
plot(f0,CA300,'r','linewidth',2);
hold off
set(gca,'yscale','log')

figure(3)
hatranK = (1)^2*400e-6*10000e-6*hatran(:,2);   %% abs coeff/cm
hatranOD = hatranK*5*1000*100;           %% OD for 5 km (assuming constant WV MR!!!
addpath /home/sergio/KCARTA/MATLAB
[d,w] = readkcstd('/home/sergio/KCARTA/UTILITY/l2s_kc120_H16.dat');
%% mimic 1x = US Std and 5x = VeryHumid
semilogy(w,sum(d'),'b',hatran(:,1),hatranOD,'g',...
         f0,1*CA300*(1)^2*400e-6*10000e-6*5*1000*100,'r',...
         f0,5*CA300*(1)^2*400e-6*10000e-6*5*1000*100,'m')
axis([2000 2600 1e-6 1e+6])	 

haha = interp1(f0,1*CA300*(1)^2*400e-6*10000e-6*5*1000*100,w);
semilogy(w,haha./sum(d'))
axis([2000 2600 1e-6 1e+6])
axis([2300 2500 1e-6 1e-1])
set(gca,'yscale','linear')
axis([2300 2500 1e-6 0.075])
title('CA co2/wv   / OD all')

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/stdNH3_1100mb_op_400ppm.rtp');
mmw = mmwater_rtp(h,p);

profile1 = load('/home/sergio/SPECTRA/IPFILES/std_gx1x_6');
profile2 = load('/home/sergio/SPECTRA/IPFILES/std_gx2x_6');
plot(profile2(:,3)./profile2(:,2));  %% 3.85e-4 ppm

T = profile1(:,4); plot(T);
stemp = T(1);

xWV  = profile1(:,3)./profile1(:,2);
xCO2 = profile2(:,3)./profile2(:,2);

rho = profile1(:,2)*101325 ./ 8.31 ./profile1(:,4);
rho = rho/1e6*6.023e23;   % molecules/cm3
rho = rho/(2.6867805e19); % amagat 

frx = [];
rad0x = [];
rad1x = [];

factor = 5

for chunk = 2355:25:2455
  addpath /home/sergio/HITRAN2UMBCLBL/FORTRAN/for2mat/
  addpath /asl/matlib/aslutil
  fname = ['/asl/data/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/CO2_385ppmv/r' num2str(chunk) '_g2.dat'];
  [gid,fr,kcomp,B] = for2mat_kcomp_reader(fname);
  odusual = B * squeeze(kcomp(:,:,6));
  odusual = odusual.^4;

  rad0 = ttorad(fr,stemp)';
  rad1 = ttorad(fr,stemp)';
  for ii = 1 : 100
    fprintf(1,'%4i %4i \n',chunk,ii);
    odx = odusual(:,ii);
    tr = exp(-odx);
    rad0 = rad0.*exp(-odx) + (1-tr).*ttorad(fr,T(ii))';

    CA = wv_co2_ca(T(ii),profile1(ii,2),fr);
    ODCA = rho(ii)*rho(ii)*xWV(ii)*xCO2(ii)*CA;
    odx = odusual(:,ii) + factor*ODCA' * 2.5e4;  %% assume path length of 250 m, this is correct at surface
                                          %% assume x factor (factor = 5 for very wet atm)  wetter than US Std
    tr = exp(-odx);
  
    rad1 = rad1.*exp(-odx) + (1-tr).*ttorad(fr,T(ii))';
  end
  frx  = [frx fr];
  rad0x=[rad0x; rad0];
  rad1x=[rad1x; rad1];  
  
  figure(1); plot(frx,rad2bt(frx,rad0x),'b',frx,rad2bt(frx,rad1x),'r',frx,ones(size(frx))*stemp,'k');
  figure(2); plot(frx,rad2bt(frx,rad1x)-rad2bt(frx,rad0x),'r')
  disp(' ')
  pause(0.1)
end

figure(1); plot(frx,rad2bt(frx,rad0x),'b',frx,rad2bt(frx,rad1x),'r',frx,ones(size(frx))*stemp,'k');
figure(2); plot(frx,rad2bt(frx,rad1x)-rad2bt(frx,rad0x),'r')

