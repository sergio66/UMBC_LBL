gid = input('Enter which gas : ');
wn = input('Enter which chunk : ');

%% for tt = 1 : 11
for tt = 6
  strx = [num2str(floor(wn),'%05d')];
  JOB = [num2str(gid,'%02d') strx num2str(tt,'%02d')];
  clust_runXtopts_savegasN_file
end 
