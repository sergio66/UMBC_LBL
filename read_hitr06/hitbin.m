
addpath '/home/motteler/radtrans/read_hitr06'

wlist = [];

for gid = 3:32

  s = read_hitran1(605, 2805, 0, gid, '/asl/data/hitran/h04.by.gas');

  wlist = [wlist; s.wnum];
  fprintf(1, 'gid=%d, length(wlist)=%d\n', gid, length(wlist));

end

