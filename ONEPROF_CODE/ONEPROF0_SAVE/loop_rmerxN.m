clear all; 

gids = input('Enter gasID list (or -1 to prompt for start/stop gasID) : ');
if gids == -1
  iA = input('Enter Start and Stop gasIDs : ');
  gids = iA(1) : iA(2);
end

for gg = 1 : length(gids)
  gid = gids(gg);

  nbox = 5;
  pointsPerChunk = 10000;
  freq_boundaries

  diroutX = [dirout '/g' num2str(gid) '.dat/'];
  diroutX = dirout;

  thedir = dir([diroutX 'std*.mat']);

  numfiles(gg) = length(thedir);
  badfiles(gg) = length(thedir);

  if length(thedir) > 0
    clear lala

    for ii = 1 : length(thedir)
      lala(ii) = thedir(ii).bytes;
      if thedir(ii).bytes <= eps
        fname = [diroutX thedir(ii).name];
        rmer = ['!/bin/rm ' fname];
        fprintf(1,'  >> rm %s \n',fname);
        eval(rmer);
      end
    end

    plot(lala);
    oo = find(lala <= eps);
    badfiles(gg) = length(oo);
    fprintf(1,'gas %2i found %5i bad files out of %5i or %6.2f chunks \n',gid,length(oo),length(thedir),length(thedir)/11)
  else
    fprintf(1,'no files in %s \n',diroutX);
  end

  pause(0.1)
end

plot(gids,numfiles,'bo-',gids,badfiles,'rx-');
hl = legend('All files','Bad files'); set(hl,'fontsize',10); grid
