cd /asl/data/hitran/H2020

fhitranver = 'h20.by.gas';
fdir0 = ['../' fhitranver '/'];
if ~exist(fdir0,'dir')
  mker = ['!mkdir -p ' fdir0];
  eval(mker)
end

fTAPE7 = '6175cfc3.par'; %% downloaded from HITRAN website

for ii = 1 : 49
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

disp('I bet g30, 35, 42 will be empty files')
disp('I bet g30, 35, 42 will be empty files')
disp('I bet g30, 35, 42 will be empty files')

iUnZip = +1;  %% if you did not bring   down the .par files
iUnZip = -1;  %% if you         brought down the .par files
if iUnZip > 0
  unzipper = ['!unzip G30_35_42/30_hit08.zip']; eval(unzipper)
  unzipper = ['!unzip G30_35_42/35_hit08.zip']; eval(unzipper)
  unzipper = ['!unzip G30_35_42/42_hit08.zip']; eval(unzipper)
end

iDoG30_35_42 = +1;
if iDoG30_35_42 > 0
  trer = ['!tr -d ''\015\032'' < G30_35_42/30_hit08.par  > ../' fhitranver '/g30.dat']; eval(trer);
  trer = ['!tr -d ''\015\032'' < G30_35_42/35_hit08.par  > ../' fhitranver '/g35.dat']; eval(trer);
  trer = ['!tr -d ''\015\032'' < G30_35_42/42_hit16.par  > ../' fhitranver '/g42.dat']; eval(trer);
end  
