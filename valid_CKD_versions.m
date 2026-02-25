%% this is for run8watercontinuum.m, do_local_lineshape_CKD.m

CKD_orig      = [0 21 23 24];
CKD_MT_orig   = [1 2 3 4 5 6];
CKD_MT_latest = [25 27 32 43];

allowedCKD = [CKD_orig CKD_MT_orig CKD_MT_latest];

if (~ismember(CKD,allowedCKD))
  disp('valid CKD versions = ')
  sort(allowedCKD)
  error('invalid CKD version! please retry!')
end
