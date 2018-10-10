function qisotopes = getq_oldVSnew(hlist_qtips,hitran_version,A,B,C,D,G,...
                                   mass_info,gasID,tempr);

%% since the latest HITRAN have updated QTIPS functions, we need to run
%% the new_q fcn which uses Lagrange interps instead of polynomial approx
%% mass_info is from load_mass_dat =>
%% mass_info = [mass_abn mass_dgn mass_QT296];

if length(intersect(hitran_version,hlist_qtips)) ~= 1
  %hitran_version
  %chlist_qtips
  if strcmp(hitran_version,'h04') == 1
    find_qnew_isotopes_H04;
  elseif strcmp(hitran_version,'h08') == 1
    find_qnew_isotopes_H08;
  elseif strcmp(hitran_version,'h12') == 1
    find_qnew_isotopes_H12;
  elseif strcmp(hitran_version,'g15') == 1
    find_qnew_isotopes_H16;
  elseif strcmp(hitran_version,'h16') == 1
    find_qnew_isotopes_H16;
  elseif strcmp(hitran_version,'h17') == 1
    %% these are fake for HITRAN LM
    find_qnew_isotopes_H16;
  elseif strcmp(hitran_version,'h18') == 1
    %% these are fake for HITRAN LM  
    find_qnew_isotopes_H16;
  else
    hitran_version
    error('looking for H04/H08/H12/H16 in new_q')
  end
  qisotopes = qfcnALL;
else
  qisotopes = [];
end

