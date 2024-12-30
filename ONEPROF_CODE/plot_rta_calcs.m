figure(1); clf
semilogy(profileT,profileP*1013); set(gca,'ydir','reverse'); ylim([0.01 1000])
  xlabel('T [K]'); ylabel('P [mb]')

[fc,qc] = quickconvolve(wall,rall,0.5,0.5);
plot(wall,rad2bt(wall,rall),fc,rad2bt(fc,qc))
if exist('sarta_calc')
  ax = axis;
  plot(wall,rad2bt(wall,rall),'b',fc,rad2bt(fc,qc),'r',h.vchan,sarta_calc,'k')
  xlim([ax(1) ax(2)])
  plot(fc,rad2bt(fc,qc),'r',h.vchan,sarta_calc,'k')
  xlim([min(fc) max(fc)])
end
xlabel('Wavenumber [cm-1]'); ylabel('BT [K]')

figure(2); clf
[fc,wgc] = quickconvolve(wall,wgtall,0.5,0.5);
pcolor(fc,1013*profileP,wgc');; shading interp; colorbar; colormap jet; %% set(gca,'ydir','reverse')
pcolor(fc,1013*profileP,wgc');; shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000])
  xlabel('Wavenumber [cm-1]'); ylabel(' P [mb]'); title('WgtX (d\tau) from RTA')

figure(3); clf
kcarta    = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/wgtfcn_jac.mat');
kcartarad = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/wgtfcn_jac.mat');
airslevels = load('/home/sergio/MATLABCODE/airslevels.dat');
airslays = plevs2plays(airslevels);
airslays = airslays(4:100);
pcolor(kcarta.fout,airslays,kcarta.jout'); colormap jet; colorbar; xlim([min(fc) max(fc)]); shading interp; set(gca,'ydir','reverse')
ylim([100 1000])
  xlabel('Wavenumber [cm-1]'); ylabel(' P [mb]'); title('WgtX (d\tau) from KCARTA')
dpkc = diff(airslays); dpkc = [dpkc; dpkc(end)*1.02]';
dpkc = ones(length(kcarta.fout),1) * -dpkc;
pcolor(kcarta.fout,airslays,kcarta.jout'./dpkc'); colormap jet; colorbar; xlim([min(fc) max(fc)]); shading interp; set(gca,'ydir','reverse')
ylim([100 1000])
  xlabel('Wavenumber [cm-1]'); ylabel(' P [mb]'); title('WgtX (d\tau)/dp from KCARTA')
 caxis([0 1]*4e-3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
plot(kcarta.fout,sum(kcarta.jout,2),'b.-',fc,sum(wgc'),'linewidth',2); xlim([wall(1) wall(end)])
  title('Sum (Weighting functionsX (sum(d\tau))'); legend('LBL','kCARTA','location','best');

figure(5); clf
[fc,contrc] = quickconvolve(wall,contrall,0.5,0.5);
pcolor(fc,1013*profileP,contrc');; shading interp; colorbar; colormap jet; %% set(gca,'ydir','reverse')
pcolor(fc,1013*profileP,contrc');; shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000])
  xlabel('Wavenumber [cm-1]'); ylabel(' P [mb]'); title('Contr functions from RTA')

figure(5); clf
dp = diff(1013*profileP); dp = [dp dp(end)*1.02];
dp = ones(length(fc),1) * dp;
wgc2 = wgc./dp;
pcolor(fc,1013*profileP,wgc2');; shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000])
  xlabel('Wavenumber [cm-1]'); ylabel(' P [mb]'); title('WgtX from RTA d\tau/dp')
  caxis([0 1]*4e-3)

figure(6); clf
iChan = [1 8]; plot(contrc(iChan,:),1013*profileP,wgc(iChan,:),1013*profileP,p.gas_1(1:100)/10e22,1013*profileP); set(gca,'ydir','reverse')
iChan = [  8]; plot(contrc(iChan,:),1013*profileP,'b.-',wgc(iChan,:),1013*profileP,'rx-',p.gas_1(1:100)/10e22,1013*profileP,'k--'); set(gca,'ydir','reverse')
iChan = [  1]; iChan = [  8]; plot(contrc(iChan,:),1013*profileP,wgc(iChan,:),1013*profileP,p.gas_1(1:100)/10e22,1013*profileP,'k--'); set(gca,'ydir','reverse')

iChan = [  1]; kcChan = find(kcarta.fout >= fc(iChan),1);
iChan = [  8]; kcChan = find(kcarta.fout >= fc(iChan),1);
plot(kcarta.jout(kcChan,:),airslays,'b.-',wgc(iChan,:),1013*profileP,'rx-'); set(gca,'ydir','reverse'); grid; ylim([700 1000]); 
  title(['Wgt fcns X for ' num2str(fc(iChan)) ' cm-1']); legend('Orig','New','location','best')
plot(kcarta.jout(kcChan,:)./dpkc(1,:),airslays,'b.-',wgc(iChan,:)./dp(1,:),1013*profileP,'rx-'); set(gca,'ydir','reverse'); grid; ylim([700 1000]); 
  title(['Wgt fcns d(\tau)/dp for ' num2str(fc(iChan)) ' cm-1']); legend('Orig','New','location','best')

figure(7); clf
tau2space = cumsum(odall,1); tau2space = exp(-tau2space);
[fc,tau2spc] = quickconvolve(wall,tau2space,0.5,0.5);
pcolor(fc,1013*profileP,tau2spc'); colormap jet; colorbar; shading interp; set(gca,'ydir','reverse')
for ii = 1 : length(fc)
  wgtnew(ii,:) = gradient(tau2spc(ii,:),-1013*profileP);
end
pcolor(fc,1013*profileP,wgtnew'); colormap jet; colorbar; shading interp; set(gca,'ydir','reverse')
caxis([0 1]*0.004); title('d\tau/dp  [1/mb]'); xlabel('Wavenumber cm-1'); ylabel('P [mb]')

figure(3); ylim([100 1000])
figure(5); ylim([100 1000])
figure(7); ylim([100 1000])

[fc,odc] = quickconvolve(wall,odall,0.5,0.5);
figure(8);
pcolor(fc,1013*profileP,log10(odc')); shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); title('log10(OD)')

[fc,radNc] = quickconvolve(wall,radNall,0.5,0.5);
figure(9);
pcolor(fc,1013*[profileP profileP(end)*1.01],(rad2bt(fc,radNc))'); shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); title('BT through the layers')

