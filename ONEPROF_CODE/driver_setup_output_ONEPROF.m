addpath /home/sergio/SPECTRA

nbox = 5;
pointsPerChunk = 10000;

glist = [1 2 3 4 5 6 8 12];

for ggx = 1 : length(glist)
  gid = glist(ggx);
  freq_boundaries;

  if ~exist(dirout)
    fprintf(1,'making %s \n',dirout)
    mker = ['!mkdir -p ' dirout];
    eval(mker)
  end
end

all_gases_REMOVE_dangerous

