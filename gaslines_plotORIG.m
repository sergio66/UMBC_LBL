t0 = 300;
wv1 = 20;   wv2 = 5000; dv = 10;
wv1 = 2850; wv2 = 3250; dv = 1;

wv = wv1:dv:wv2;
tv = ttorad(wv,t0);

gasids = ['H2OCO2 O3N2O COCH4'];
gasname = ['H2OCO2xO3N2OxCOCH4'];
for gg = 1 : 26
  gasID = gg;

  %%fnamePRE='/asl/data/hitran/h92.by.gas/g';
  %fnamePRE='/salsify/scratch4/h96.by.gas/g';

  fnamePRE='/asl/data/hitran/h98.by.gas/g';
  fnamePRE='/asl/data/hitran/h92.by.gas/g';
  do_HITRAN_vers;
  fnamePRE = [HITRAN '/g'];; 
  
  fnamePOST='.dat';
  fnameIN=int2str(gasID);
  hitlin_fname=[fnamePRE fnameIN fnamePOST];

  start = wv1;
  stop  = wv2;

  line=hitread(start,stop,0,gasID,hitlin_fname);

%  ii = find(line.wnum >= wv1 & line.wnum <= wv2);
%  semilogy(line.wnum,line.stren,'.',line.wnum(ii),line.stren(ii),'r.',...
%           wv,tv/max(tv)*max(line.stren));

  if line.linct > 0
    clf
    semilogy(line.wnum,line.stren,'.')
    %jj = (gg-1)*3 + (1:3);
    %title(['gas = ' gasids(jj)]);
    title(['gas = ' num2str(gg)]);
    pause
    %fname = ['gasid' gasname(jj)];
    %printfig(1,fname,'jpg');
    end
  end

%plotyy(line.stren,line.wnum,line.stren(ii),10000./line.wnum(ii),'semilogx'); 
%plotyy(line.stren,line.wnum,line.stren(ii),line.wnum(ii),'semilogx'); 
