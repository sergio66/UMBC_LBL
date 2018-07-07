function  qfcn = new_q(tempr,mass_info,gasID,line,hitran_version);

%%% this is not so fast as it does a lot of irrelevant but needed I/O
%%% slows you down a lot when you loop over near,medium,far meshes

%% since the latest HITRAN have updated QTIPS functions, we need to run
%% the new_q fcn which uses Lagrange interps instead of polynomial approx
%% mass_info is from load_mass_dat =>
%% mass_info = [mass_abn mass_dgn mass_QT296];

if strcmp(hitran_version,'h04') == 1
  find_qnew_isotopes_H04;
elseif strcmp(hitran_version,'h08') == 1
  %find_qnew_isotopes_H08_OLD;
  find_qnew_isotopes_H08;
else
  hitran_version
  error('looking for H04 or H08 in new_q')
  end

for ii = 1 : x
  yoyo = find(line.iso == ii);
  qfcn(yoyo) = qfcnALL(ii);
  end
qfcn = qfcn';


