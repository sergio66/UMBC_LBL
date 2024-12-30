addpath /home/sergio/SPECTRA
addpath /asl/matlib//h4tools
addpath /asl/matlib//aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

set_file_names

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

if isfield(p,'rcalc')
  sarta_calc = rad2bt(h.vchan,p.rcalc);
end

printarray([junkone.rtpProf.mpres*1013  pxP junkP'*1013   junkone.rtpProf.mtemp pxT junkT'])
fprintf(1,'read in data from %s and %s \n',frtp,dirSAVE)

junk = input('Enter +1/default to proceed, -1 to stop : ');
if length(junk) == 0
  junk = +1;
end

if junk < 0
  error('stopping')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stemp = 304.0

nbox = 5;
pointsPerChunk = 10000;
gid = 1; freq_boundaries

wall     = [];
rall     = [];
radNall  = [];
wgtall   = [];
contrall = [];
odall    = [];

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
  tsurf = stemp;

  rad = ttorad(w,tsurf);
  radN(length(profileT)+1,:) = rad;
  for ll = length(profileT) : -1 : 1
    T = profileT(ll);
    tau = od(ll,:);
    x = exp(-tau);
    rad = rad .* x + ttorad(w,T).*(1-x);
    radN(ll,:) = rad;
    if ll > 1
      wah = od(1:ll-1,:);
      wah = sum(wah,1);
      wgt(ll,:) = (1-x).*exp(-wah);

      contr(ll,:) = ttorad(w,T).*(1-x) .*exp(-wah);

      %figure(1); plot(w,exp(-wah)); ylabel('T(ll->TOA'); title(num2str(ll)); 
      %figure(2); plot(w,wgt(ll,:)); ylabel('W G T ');    title(num2str(ll)); 
      %pause(0.1)
    else
      wgt(ll,:) = (1-x);
      contr(ll,:) = ttorad(w,T).*(1-x);
    end
  end
  contr = contr./(ones(length(profileT),1)*rad);
  wall     = [wall    w];
  rall     = [rall    rad];
  radNall  = [radNall radN];
  wgtall   = [wgtall   wgt];
  contrall = [contrall contr];
  odall    = [odall    od];

  figure(1); clf; plot(w,rad2bt(w,rad))
  figure(2); clf; pcolor(w,profileP,wgt); set(gca,'ydir','reverse'); shading interp; colorbar

  %disp('ret to continue'); pause;
  pause(0.1)

end       %% loop over ff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_rta_calcs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSave = input('save (-1/+1) : ');
if length(iSave) == 0
  iSave = -1;
end
if iSave > 0
  if findstr(dirSAVE,'ChrisLayers')
    saver = ['save ' dirSAVE '/pbl100_simple_rta.mat wall rall fc odc radNc pxP pxT'];
  elseif findstr(dirSAVE,'AIRSLayers')
    saver = ['save ' dirSAVE '/airs100_simple_rta.mat wall rall fc odc radNc pxP pxT'];
  end
  eval(saver)
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSave = input('show summary of PBL and AIRS rad transfer? (-1/+1) : ');
if iSave > 0
  airsrad = load('/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/SAVE/PBLTEST/AIRSLayers/airs100_simple_rta.mat');
  pblrad  = load('/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/SAVE/PBLTEST/ChrisLayers/pbl100_simple_rta.mat');
  
  figure(1); clf; plot(airsrad.wall,rad2bt(airsrad.wall,airsrad.rall),airsrad.wall,rad2bt(airsrad.wall,pblrad.rall))
    title('BT TOA (b) AIRS (r) PBL')
  figure(2); clf; plot(airsrad.wall,rad2bt(airsrad.wall,airsrad.rall)-rad2bt(airsrad.wall,pblrad.rall))
    title('BTD TOA AIRS-PBL')
  
  figure(3); clf; pcolor(airsrad.fc,1:100,log10(pblrad.odc'./airsrad.odc')); ylabel('log10(PBL OD/AIRS OD)'); 
    shading interp; title('Expect different ODs \newline since different layering'); colorbar
    set(gca,'ydir','reverse')
  figure(4); clf; pcolor(airsrad.fc,1:101,rad2bt(airsrad.fc,airsrad.radNc)'-rad2bt(airsrad.fc,pblrad.radNc)'); ylabel('BTD AIRS-PBL'); 
    shading interp; title('Expect different BTs in middle atm \newline since different layering'); colorbar
    set(gca,'ydir','reverse')

  figure(5); clf; semilogy(airsrad.pxT,airsrad.pxP,pblrad.pxT,pblrad.pxP); set(gca,'ydir','reverse'); xlabel('T [K]'); ylabel('P [mb]'); ylim([0 1000])
  figure(6); clf; semilogy(airsrad.pxT,1:100,pblrad.pxT,1:100); set(gca,'ydir','reverse'); xlabel('T [K]'); ylabel('layer number')
  figure(7); clf; semilogy(airsrad.pxT-pblrad.pxT,1:100); set(gca,'ydir','reverse'); xlabel('AIRS-PBL \delta T [K]'); ylabel('layer number'); plotaxis2;
end
