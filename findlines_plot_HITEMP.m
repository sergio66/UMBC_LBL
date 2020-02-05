function [iYes,line] = findlines_plot_HITEMP(wv1,wv2,gid,HITRAN)

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
% tv= ttorad(wv,t0);

line = [];

gasID = gid;
hstr = num2str(HITRAN-2000,'%02d');

fnamePRE = ['/asl/data/hitran/HITEMP/h' hstr '.by.gas/g'];

fnamePOST    = '.dat';
fnameIN      = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];
fprintf(1,'in findlines_plot.m we are opening %s for lineparams for gasID %2i \n',hitlin_fname,gasID)

start = wv1;
stop  = wv2;

iYes = -1;

clf
if gid ~= 7 & gid ~= 22
  %% usual gases
  line = hitread(start,stop,0,gasID,hitlin_fname,-1);

  if line.linct > 0
    iYes = +1;
    semilogy(line.wnum,line.stren,'.'); title(['gasID = ' num2str(gasID)])
  end

elseif gid == 7
  %% O2
  line = hitread(start,stop,0,gasID,hitlin_fname,-1);
  if line.linct > 0
    iYes = +1;
    semilogy(line.wnum,line.stren,'.'); title(['gasID = ' num2str(gasID)])
  end
  %% also add on continuum check as needed
  if wv1 <= 1850 & wv2 >= 1340
    iYes = 1;
  end

elseif gid == 22
  line=hitread(start,stop,0,gasID,hitlin_fname,-1);
  %% N2
  if line.linct > 0
    iYes = +1;
    semilogy(line.wnum,line.stren,'.'); title(['gasID = ' num2str(gasID)])
  end
  %% also add on continuum check as needed
  if wv1 <= 350 & wv2 >= 0
    iYes = 1;
  end
  %% also add on continuum check as needed
  if wv1 <= 2670 & wv2 >= 2085
    iYes = 1;
  end


end

