function [iYes,line] = findlines_plot(wv1,wv2,gid,HITRANyear,HorG);

if nargin < 4
  HITRANyear = 2012;
  HITRANyear = 2016;
  HITRANyear = 2020;
  HorG = +1;      %% HITRANyear or GEISA = 1 ==> HITRANyear        HITRANyear or GEISA = -1 ==> GEISA
elseif nargin < 5
  HorG = +1;      %% HITRANyear or GEISA = 1 ==> HITRANyear        HITRANyear or GEISA = -1 ==> GEISA
end

if abs(HorG) ~= 1
  error('HorG = +/-1 for HITRAN or GEISA')
end

if HorG == 1
  if HITRANyear == 1992
    hstr = '92';
  elseif HITRANyear == 1996
    hstr = '96';
  elseif HITRANyear == 1998
    hstr = '98';
  elseif HITRANyear == 2000
    hstr = '2k';
  elseif HITRANyear > 2000
    hstr = num2str(HITRANyear-2000,'%02d');
  elseif HITRANyear == 2004
    hstr = '04';    
  elseif HITRANyear == 2008
    hstr = '08';    
  elseif HITRANyear == 2012
    hstr = '12';    
  elseif HITRANyear == 2016
    hstr = '16';    
  elseif HITRANyear == 2020
    hstr = '20';    
  elseif HITRANyear == 2024
    hstr = '24';
  else
    HITRANyear
    error('only have HITRAN 96,98,2k,04, ... 24')
  end
elseif HorG == -1
  if HITRANyear == 2015
    hstr = '15';
  else
    HITRANyear
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
  fnamePRE = ['/umbc/xfs3/strow/asl/rta/hitran/h' hstr '.by.gas/g'];  
  fnamePRE = [do_HITRAN_vers '/g'];
  fnamePRE = [hitranpath '/h' hstr '.by.gas/g'];    
  
else
  fnamePRE = ['/asl/data/geisa/g' hstr '.by.gas/g'];
  fnamePRE = ['/asl/rta/geisa/g' hstr '.by.gas/g'];
  fnamePRE = ['/umbc/xfs3/strow/asl/rta/geisa/h' hstr '.by.gas/g'];
  fnamePRE = [do_GEISA_vers '/g'];
  fnamePRE = [geisapath '/g' hstr '.by.gas/g'];      
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

