function chunklen = find_chunklen(fstep,nbox,ffin);

% function chunklen = find_chunklen(fstep,nbox,ffin);
chunklen0 = fstep/(nbox*ffin);
chunklenF = floor(fstep/(nbox*ffin));
chunklenC = ceil(fstep/(nbox*ffin));
chunklenX = fix(fstep/(nbox*ffin));

if ((abs(chunklen0 - chunklenF) < 0.1) & (abs(chunklen0 - chunklenX) < 0.1))
  chunklen = chunklenX;
else
  disp('chunklen0,chunklenF,chunklenC,chunklenX = ')
  chunklen0
  chunklenF
  %chunklenC
  chunklenX
  error('error find_chunklen');
  end
