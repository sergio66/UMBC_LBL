function qfcn = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                           qisotopes,mass_info,gasID,line,tempr);

if length(intersect(hitran_version,hlist_qtips)) == 1
  %%% this is fast
  qfcn = q(A,B,C,D,G,line,tempr);
else
  %% this is very fast as qisotopes has already been calculated for this layer
  [x,y] = size(mass_info);
  %%% for safety sake, set qfcn = 1 for ALL lines
  %%% this helps gas27, which has 2 isotopes in H08 version, but only
  %%%   the main isotope has partition fcn info 
  %%%   (qtips : see Global_Data_HITRAN2008/molparam.txt)
  qfcn = ones(size(line.iso));
  for ii = 1 : x
    yoyo = find(line.iso == ii);
    qfcn(yoyo) = qisotopes(ii);
    end
  qfcn = qfcn';
  end

