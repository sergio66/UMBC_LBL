
% compare matlab and mexized C HITRAN readers

more off

hdat = '/asl/data/hitran/h04.by.gas';
v1 = 600;
v2 = 2800;

for gid = 1 : 39

  fprintf(1, 'reading data for gas %d...\n', gid)
  tic; L1 = read_hitran2(v1, v2, 0, gid, hdat); toc
  tic; L2 = read_hitran(v1, v2, 0, gid, hdat); toc

  % fprintf(1, 'comparing data for gas %d...\n', gid)

  if isempty(L1.igas) & isempty(L2.igas)
    continue
  end

  if ~isequal(L1, L2)
    fprintf(1, '*** WARNING readers disagree for gas %d ***\n', gid)
  end

end

