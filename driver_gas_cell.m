function [w,d] = driver_gas_cell(gasID)

cd /home/sergio/git/SPECTRA/

%% gasID = input('Enter gasID : ');

if gasID == 1
  tempfilename = '/home/sergio/git/SPECTRA/IPFILES/xgas_cell_water';   %%% edit this for your path/file name
  gas_cell_water
  d = d8conkc+d8kc;
%elseif gasID == 2
%  tempfilename = '/home/sergio/git/SPECTRA/IPFILES/xgas_cell_co2';    %%% edit this for your path/file name
%  gas_cell_co2
%  d = d8;
elseif gasID > 1
  if gasID == 2
    disp('WARNING : should really be running gas_cell_co2')
  end
  tempfilename = '/home/sergio/git/SPECTRA/IPFILES/xgas_cell';         %%% edit this for your path/file name
  gas_cell_others
end
