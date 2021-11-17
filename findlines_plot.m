function [iYes,line] = findlines_plot(wv1,wv2,gid,HITRAN,HorG);

if nargin < 4
  HITRAN = 2012;
  HITRAN = 2016;
  HITRAN = 2020;
  HorG = +1;      %% HITRAN or GEISA = 1 ==> HITRAN        HITRAN or GEISA = -1 ==> GEISA
elseif nargin < 5
  HorG = +1;      %% HITRAN or GEISA = 1 ==> HITRAN        HITRAN or GEISA = -1 ==> GEISA
end

if abs(HorG) ~= 1
  error('HorG = +/-1 for HITRAN or GEISA')
end

if HorG == 1
  if HITRAN == 1992
    hstr = '92';
  elseif HITRAN == 1996
    hstr = '96';
  elseif HITRAN == 1998
    hstr = '98';
  elseif HITRAN == 2000
    hstr = '2k';
  elseif HITRAN > 2000
    hstr = num2str(HITRAN-2000,'%02d');
  end
elseif HorG == -1
  if HITRAN == 2015
    hstr = '15';
  else
    error('only have GEISA 2015')
  end
end

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
% tv= ttorad(wv,t0);

line = [];

gasID = gid;

%%fnamePRE = '/asl/data/hitran/h92.by.gas/g';
%fnamePRE = '/salsify/scratch4/h96.by.gas/g';

if HorG == 1
  fnamePRE = ['/asl/data/hitran/h' hstr '.by.gas/g'];
  fnamePRE = ['/asl/rta/hitran/h' hstr '.by.gas/g'];
else
  fnamePRE = ['/asl/data/geisa/g' hstr '.by.gas/g'];
  fnamePRE = ['/asl/rta/geisa/g' hstr '.by.gas/g'];
end

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

