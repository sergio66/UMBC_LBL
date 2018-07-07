function [outwave,out_array] = driver_run8_parallel(gas,wn1,wn2,profname,dv,options);


matlabpool open
poolSize = matlabpool('size');
if poolSize == 0
   error('parallel:demo:poolClosed .... This demo needs an open MATLAB pool to run.');
end
fprintf('This demo is running on %d MATLABPOOL workers.\n', matlabpool('size'));

%% assumes dv cm-1 chunks = 10000 pts
boo = wn1 : dv : wn2

outwave = [];
out_array = [];
parfor hh = 1 : length(boo)-1
  index = (1:10000) + (hh-1)*10000;
  wn = boo(hh);
  if nargin == 6
    [g(hh) wx odx] =  run8_parallel(gas,boo(hh),boo(hh)+dv,profname,options);
  else
    [g(hh) wx odx] =  run8_parallel(gas,boo(hh),boo(hh)+dv,profname);
  end
  saver = ['save /home/sergio/HITRAN2UMBCLBL/MARS_MAKEIR/Develop/mars_co2_ods_' num2str(wn) '.mat wx odx'];
  eval(saver)
end

matlabpool close

