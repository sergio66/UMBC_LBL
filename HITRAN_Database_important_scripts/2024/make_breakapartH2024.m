%cd /asl/data/hitran/H2024
cd /umbc/xfs3/strow/asl/rta/hitran/H2024

fhitranver = 'h24.by.gas';
fdir0 = ['../' fhitranver '/'];
if ~exist(fdir0,'dir')
  mker = ['!mkdir -p ' fdir0];
  eval(mker)
end

fTAPE7 = '6175cfc3.par'; %% H2020 downloaded from HITRAN website
fTAPE7 = '6984229e.par'; %% H2024 downloaded from HITRAN website, Feb 5, 2026

for ii = 1 : 61
  fprintf(1,' making file %2i \n',ii)
  fOUT = [fdir0 '/g' num2str(ii) '.dat'];  
  
  if ii < 10
    grepper = ['!grep "^ ' num2str(ii) '[1-9]" ' fTAPE7 ' > junk.lines'];
  else
    grepper = ['!grep "^'  num2str(ii) '[1-9]" ' fTAPE7 ' > junk.lines'];
  end
  eval(grepper)

  %% remove junk chars 
  trer = ['!tr -d ''\015\032'' < junk.lines  > ' fOUT];
  eval(trer);
end

disp('I bet g30, 35, 42, 55 will be empty files')
disp('I bet g30, 35, 42, 55 will be empty files')
disp('I bet g30, 35, 42, 55 will be empty files')

error('download  G30,35,42,55 separately to desktop, scp to chip, then run the below')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iUnZip = -1;  %% if you         brought down the .par files
iUnZip = +1;  %% if you did not bring   down the .par files
if iUnZip > 0
  unzipper = ['!unzip 30_hit24.zip']; eval(unzipper)
  unzipper = ['!unzip 35_hit24.zip']; eval(unzipper)
  unzipper = ['!unzip 42_hit24.zip']; eval(unzipper)
  unzipper = ['!unzip 55_hit24.zip']; eval(unzipper)  
end

iDoG30_35_42 = +1;
if iDoG30_35_42 > 0
  trer = ['!tr -d ''\015\032'' < 30_hit24.par  > ../' fhitranver '/g30.dat']; eval(trer);
  trer = ['!tr -d ''\015\032'' < 35_hit24.par  > ../' fhitranver '/g35.dat']; eval(trer);
  trer = ['!tr -d ''\015\032'' < 42_hit24.par  > ../' fhitranver '/g42.dat']; eval(trer);
  trer = ['!tr -d ''\015\032'' < 55_hit24.par  > ../' fhitranver '/g55.dat']; eval(trer);  
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
