%%%%%%%%%%%%%%%%%%%%%%%%%  y1ser.m %%%%%%%%%%%%%%%%%%%%%%%%
%
%   This program simply calls y1s.m with the proper variables to calculate
%   the first order mixing coefficients.
%

disp('computing 1st order mixing coefficients')
for i=1:length(pick_no)
ni=num2str(i);
eval(['ymixr(:,' ni ')=y1s(jr,freqr,elowerr,strenr,W_plus_r' ni ',''R'')''; '])
eval(['ymixp(:,' ni ')=y1s(jp,freqp,elowerp,strenp,W_plus_p' ni ',''P'')''; '])
end


