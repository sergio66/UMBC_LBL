addpath /home/sergio/SPECTRA

%% addpath '/home/motteler/radtrans/read_hitr06'

wlist = [];

for gid = 3:32

  %s1 = read_hitran1(605, 2805, 0, gid, '/asl/data/hitran/h04.by.gas');
  %wlist = [wlist; s1.wnum];
  
  s = read_hitran(605, 2805, 0, gid, [hitranpath '/h24.by.gas/g' num2str(gid) '.dat']);
  wlist = [wlist; s1.wnum];
  
  fprintf(1, 'gid=%d, length(wlist)=%d \n', gid, length(wlist));

end

