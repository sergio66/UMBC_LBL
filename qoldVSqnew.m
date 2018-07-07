function qfcn = qoldVSqnew(hitran_version,A,B,C,D,G,...
                           mass_info,gasID,line,tempr);

if length(intersect(hitran_version,{'h92','h96','h98','h2k'})) == 1
  %%% this is fast
  qfcn = q(A,B,C,D,G,line,tempr);
else
  %%% this is not so fast as it does a lot of irrelevant but needed I/O
  %%% slows you down a lot when you loop over near,medium,far meshes
  %%% only called by findUnionNew.m
  qfcn = new_q(tempr,mass_info,gasID,line,hitran_version);
  end
