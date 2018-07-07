disp('see /home/sergio/KCARTA/UTILITYBOX for : ')
disp(' (a) watercontinuum_nir_vis.m shows you the water continuum')
disp(' (b) umbclbl_2_kcarta_nir_vis.m takes this, polus continuum and makes'
disp('     f77 binary files that kcarta can read')

nbox = 5;
pointsPerChunk = 10000;
gases = [1 2 3 4 5 6 7];

%%%% to do things for MODIS vis
wnVIS = [4650 6095 8065 11550 15150 18150 20055 21250];
wnVIS = [4650 6095 8065 11550 15150 18150 20055 21250] - 12.5;
wnVIS1 = 4500;
wnVIS1 = 6050;
wnVIS2 = 22000;

%%% sims for Jack Kumer
wnVIS1 = 4000;
wnVIS2 = 4500;

fmin = wnVIS1; 
topts = runXtopts_params_smart(fmin); 
dv = topts.ffin*nbox*pointsPerChunk;
while fmin <= wnVIS2
  topts = runXtopts_params_smart(fmin);
  dv = topts.ffin*nbox*pointsPerChunk;
  fmax = fmin + dv;
  wvis0 = [];
  dvis0 = zeros(100,10000);
  for gg = 1 : length(gases)
    gasid = gg;  
    fip = ['IPFILES/trp_g' num2str(gg)];

    if gasid == 1
      [w,d] = run8water(gasid,fmin,fmax,fip,topts);  
    else
      [w,d] = run8(gasid,fmin,fmax,fip,topts);  
      end

    wvis0 = w;
    dvis0 = dvis0 + d;
    fout = ['/carrot/s1/sergio/RUN8_VISDATABASE/trp' num2str(fmin) '.mat'];
    saver = ['save ' fout ' wvis0 dvis0 '];
    eval(saver);
    end
  fout = ['/carrot/s1/sergio/RUN8_VISDATABASE/trp' num2str(fmin) '.mat'];
  saver = ['save ' fout ' wvis0 dvis0 '];
  eval(saver);
  fmin = fmin + dv;
  end

