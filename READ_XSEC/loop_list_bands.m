
% step thru xsec gas IDS and tabulate wavenumber bands

% loop on gas IDs
for gid = 51:81

  xs = read_xsec(gid);

  [nrec, nband] = size(xs);
  
  for b = 1:nband
  
    v1 = max([xs(:,b).v1]);
    v2 = min([xs(:,b).v2]);
  
    fprintf(1, '%4d %4d %10.3f %10.3f\n', gid, b, v1, v2);

  end % band loop

end % gas loop

