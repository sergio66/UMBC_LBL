addpath /home/sergio/SPECTRA

frtp    = '/home/chepplew/projects/klayers_wrk/regr49_pbl.op.rtp'; iProf = 1;
junkout = '/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/g1.dat/profH2O1205_1_6_2.mat';
junkout = '/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/g2.dat/prof1230_2_6.mat';
fone    = 'PROFILES/oneprof_regr49_pbl_1.mat';

dirSAVE = '/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/SAVE/PBLTEST/ChrisLayers/';
glist = [1 2 3 4 5 6 9 12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junkone = load(fone);

[h,ha,p,pa] = rtpread(frtp);
[h,p] = subset_rtp_allcloudfields(h,p,[],[],iProf);
p = make_rtp_plays(p);
pxP = p.plays(1:p.nlevs-1);
pxT = p.ptemp(1:p.nlevs-1);

junk = load(junkout);
junkP = junk.profile(2,:);
junkT = junk.profile(4,:);

printarray([junkone.rtpProf.mpres*1013  pxP junkP'*1013   junkone.rtpProf.mtemp pxT junkT'])

junk = input('Enter +1/default to proceed, -1 to stop : ');
if length(junk) == 0
  junk = +1;
end

if junk < 0
  error('stopping')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbox = 5;
pointsPerChunk = 10000;
gid = 1; freq_boundaries

wall   = [];
rall   = [];
wgtall = [];


iF = 0;
for freqchunk = wn1 : 25 : wn2
  iF = iF + 1;
  od = [];
  wgt = [];

  for gg = 1 : length(glist)
    gid = glist(gg);

    if gid == 1
      filefound(iF,gg) = +1;
      fin = [dirSAVE '/profH2O' num2str(freqchunk) '_1_6_2.mat'];
      x = load(fin);
      profileP = x.profile(2,:);
      profileT = x.profile(4,:);
      od = x.d;
      w  = x.w;

     toptsC.divide = -1; 

     iSortPtotal_Hi2Lo = +1;
     iSortPtotal_Hi2Lo = -1;
     [outwave,out_array] = run8watercontinuum_rtp(1,freqchunk,freqchunk+25,frtp,toptsC,iProf,iSortPtotal_Hi2Lo);

     ichan = find(sum(od,1) == max(sum(od,1)));
     semilogx(od(:,ichan),1:100,out_array(:,ichan)*100,1:100); legend('WV','WV continuum','location','best'); set(gca,'ydir','reverse'); 
     loglog(od(:,ichan),profileP*1013,out_array(:,ichan)*100,profileP*1013); legend('WV','WV continuum','location','best'); set(gca,'ydir','reverse'); 
       pause(0.1)

     od = od + out_array;

    else
      fin = [dirSAVE '/prof' num2str(freqchunk) '_' num2str(gid) '_6.mat'];
      if exist(fin)
        filefound(iF,gg) = +1;
        x = load(fin);
        profileP = x.profile(2,:);
        profileT = x.profile(4,:);
        od = od + x.d;
      else
        filefound(iF,gg) = 0;
      end

    end   %% gid = 1 or 2,3,4,5 ...
  end     %% loop over gg

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tsurf = profileT(end);
  rad = ttorad(w,tsurf);
  for ll = length(profileT) : -1 : 1
    T = profileT(ll);
    tau = od(ll,:);
    x = exp(-tau);
    rad = rad .* x + ttorad(w,T).*(1-x);
    if ll > 1
      wah = od(1:ll-1,:);
      wah = sum(wah,1);
      wgt(ll,:) = (1-x).*exp(-wah);
      %figure(1); plot(w,exp(-wah)); ylabel('T(ll->TOA'); title(num2str(ll)); 
      %figure(2); plot(w,wgt(ll,:)); ylabel('W G T ');    title(num2str(ll)); 
      %pause(0.1)
    else
      wgt(ll,:) = (1-x);
    end
  end
  wall = [wall w];
  rall = [rall rad];
  wgtall = [wgtall wgt];
  figure(1); clf; plot(w,rad2bt(w,rad))
  figure(2); clf; pcolor(w,profileP,wgt); set(gca,'ydir','reverse'); shading interp; colorbar

  %disp('ret to continue'); pause;
  pause(0.1)

end       %% loop over ff

figure(1); clf
semilogy(profileT,profileP*1013); set(gca,'ydir','reverse'); ylim([0.01 1000])

[fc,qc] = quickconvolve(wall,rall,0.5,0.5);
plot(wall,rad2bt(wall,rall),fc,rad2bt(fc,qc))

figure(2); clf
[fc,wgc] = quickconvolve(wall,wgtall,0.5,0.5);
pcolor(fc,1013*profileP,wgc');; shading interp; colorbar; colormap jet; %% set(gca,'ydir','reverse')
pcolor(fc,1013*profileP,wgc');; shading interp; colorbar; colormap jet; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000])

figure(3); clf
kcarta    = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/wgtfcn_jac.mat');
kcartarad = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/wgtfcn_jac.mat');
airslevels = load('/home/sergio/MATLABCODE/airslevels.dat');
airslays = plevs2plays(airslevels);
airslays = airslays(4:100);
pcolor(kcarta.fout,airslays,kcarta.jout'); colormap jet; colorbar; xlim([min(fc) max(fc)]); shading interp; set(gca,'ydir','reverse')
ylim([100 1000])
