function [] = compare_linemix_coeffs(dir1,dir2);

[dir1 '/*.mat']
[dir2 '/*.mat']

thedir1 = dir([dir1 '/*.dat']);
thedir2 = dir([dir2 '/*.dat']);

fprintf(1,'length(dir1) = %3i length(dir2) = %3i \n',length(thedir1),length(thedir2))

if length(dir2) < length(dir1)
  xdir = thedir2;
else
  xdir = thedir1;
end

iMiss = 0;
for ii = 1 : length(xdir)
  fname1 = [dir1 '/' xdir(ii).name];
  fname2 = [dir2 '/' xdir(ii).name];
  if ~exist(fname1)
    iMiss = iMiss + 1;
    fprintf(1,'%3i %3i %s DNE \n',ii,iMiss,fname1);
  elseif ~exist(fname2)
    iMiss = iMiss + 1;  
    fprintf(1,'%3i %3i %s DNE \n',ii,iMiss,fname2);
  else
    %wow1 = load(fname1)
    %wow2 = load(fname2)
    differ = ['!diff ' fname1 ' ' fname2];
    eval(differ)
  end
end

    