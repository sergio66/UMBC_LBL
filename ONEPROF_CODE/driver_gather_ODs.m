nbox = 5;
pointsPerChunk = 10000;

dirSAVE = '/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/SAVE/PBLTEST/ChrisLayers/';

lser = ['ls -lt ' dirSAVE ' | wc -l'];
[status,cmd] = system(lser);
cmd = str2num(cmd);
if cmd > 0
  lser = ['!ls -lt ' dirSAVE];
  eval(lser)
  disp('>>>>>')
  fprintf(1,'there are %5i files in %s \n',cmd,dirSAVE)
  iY = input('the outdir already has files, proceed (-1/DEFAULT +1) : ');
  if length(iY) == 0
    iY = -1;
  end
  if iY < 0
    error('stopping')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glist = [1 2 3 4 5 6 9 12];
for gg = 1 : length(glist)
  gid = glist(gg); freq_boundaries
  thedir = dir([dirout '/*.mat']);
  if length(thedir) > 0
    cper = ['!/bin/cp -a ' dirout '/*.mat ' dirSAVE '/.'];
    eval(cper)
  end
end

