disp(' WARNING : this will ERASE all .mat files you have made!!!!')
disp(' WARNING : this will ERASE all .mat files you have made!!!!')
disp(' WARNING : this will ERASE all .mat files you have made!!!!')
disp(' WARNING : this will ERASE all .mat files you have made!!!!')

nbox = 5;
pointsPerChunk = 10000;
gid = 1;
  freq_boundaries

fprintf(1,'Will be removing files from eg %s \n',dirout);
iYes = input('Proceed anyway??? (-1/+1) ??? ');
if iYes < 0
  return
end

fprintf(1,'Will be removing files from eg %s \n',dirout);
iYes = input('Proceed anyway??? Yes (-1/+1) ??? ');
if iYes < 0
  return
end

fprintf(1,'Will be removing files from eg %s \n',dirout);
iYes = input('Proceed anyway??? Really (-1/+1) ??? ');
if iYes < 0
  return
end

fprintf(1,'Will be removing files from eg %s \n',dirout);
iYes = input('Proceed anyway??? Are you sure (-1/+1) ??? ');
if iYes < 0
  return
end

%% empty out individual dirs
for gid = 1 : 110
  freq_boundaries
  thedir = dir([dirout 'std*.mat']);
  if length(thedir) > 0
    fprintf(1,'emptying out dir for gasID %3i \n',gid);
    rmer = ['!/bin/rm ' dirout 'std*.mat'];
    eval(rmer);
  end
end

