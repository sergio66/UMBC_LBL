function [iYes,line] = findlines_plot(wv1,wv2,gid,HITRAN);

if nargin < 4
  HITRAN = 2012;
  HITRAN = 2016;
end

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

t0 = 300;
dv = 10;
dv = 1;

wv = wv1:dv:wv2;
% tv= ttorad(wv,t0);

line = [];

gasID = gid;

%%fnamePRE = '/asl/data/hitran/h92.by.gas/g';
%fnamePRE = '/salsify/scratch4/h96.by.gas/g';

fnamePRE = ['/asl/data/hitran/h' hstr '.by.gas/g'];

fnamePOST    = '.dat';
fnameIN      = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];

start = wv1;
stop  = wv2;

iYes = -1;

clf
if gid ~= 7 & gid ~= 22
  %% usual gases
  line = hitread(start,stop,0,gasID,hitlin_fname,-1);

  if line.linct > 0
    iYes = +1;
    semilogy(line.wnum,line.stren,'.')
  end

elseif gid == 7
  %% O2
  line = hitread(start,stop,0,gasID,hitlin_fname,-1);
  if line.linct > 0
    iYes = +1;
    semilogy(line.wnum,line.stren,'.')
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
    semilogy(line.wnum,line.stren,'.')
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

