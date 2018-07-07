function [hh,ll] = band_end_find(gasID,homepath,thebands,vers,strengthM,fname,hitran_version);
%this function finds the start/stop points of the bands, so that they can
%be included as bands in the computation

% first see if file exists with ALL bands!!!!!!!!!!!!!!
% fid = fopen('CO2_COMMON/endpts.dat','r');

if gasID == 2
  filename1 = [homepath 'CO2_COMMON/endpts.dat'];
elseif gasID == 5
  filename1 = [homepath 'CO_COMMON/endpts.dat'];
end

fprintf(1,'in band_end_find.m, trying to find %s \n',filename1);
thedir = dir(filename1);
if length(thedir) == 1
  thedir
else
  disp('could not find the file ... fid == -1')
end

fid = fopen(filename1,'r');
ok = -1;

if fid ~= -1
  ab = fread(fid,inf,'float64');
  fclose(fid);
  ok = 1;

  %first read off the strength that this file was made for
  aa = ab(1:3)';       
  strn = aa(1);
  if strn <= strengthM          %things seem to be ok
    len = length(ab); 
    bb = ab(4:len); 
    bb = reshape(bb,length(bb)/3,3);

    file_bands = bb(:,1);
    ll = bb(:,2)';
    hh = bb(:,3)';

    %now check that all of the_bands are in file_bands;
    band_diff = setdiff(thebands,file_bands);
    if length(band_diff > 0) 
      ok = -1;                   %file does NOT have all bands
    end      

    %now see if there are TOO many bands in the data file!!!!!
    band_diff = setdiff(file_bands,thebands);
    if length(band_diff > 0)    %too many bands in the file!!!!!!    
      ism = ismember(file_bands,thebands);
      ism = find(ism);
      ll = ll(ism);
      hh = hh(ism);
    end

  else                          %file is for larger stren, so we might miss
    ok = -1;                    %so redo things!!!!!!!!!
  end
end

if ok < 0
  fprintf(1,'data file not found or incomplete -- have to create one!! \n');
  fprintf(1,' might take a few seconds ...... \n');

  [hh,ll] = band_end(gasID,thebands,vers,strengthM,fname,hitran_version);

  if gasID == 2
    filename1 = [homepath 'CO2_COMMON/endpts.dat'];
    lala = exist([homepath '/CO2_COMMON'],'dir');
  elseif gasID == 5
    filename1 = [homepath 'CO_COMMON/endpts.dat'];
    lala = exist([homepath '/CO_COMMON'],'dir');
  end
  
  if lala == 0
    if gasID == 2
      disp(' >> band_end_find.m is basically looking for /home/sergio/SPECTRA/CO2_COMMON/endpts.dat');
      disp(' >> and either cannot find the file or write to the dir');
      disp(' >> chances are you are NOT in /home/sergio/SPECTRA');
      fprintf(1,'%s \n',[homepath 'CO2_COMMON/endpts.dat'])      
    elseif gasID == 5
      disp(' >> band_end_find.m is basically looking for /home/sergio/SPECTRA/CO2_COMMON/endpts.dat');
      disp(' >> and either cannot find the file or write to the dir');
      disp(' >> chances are you are NOT in /home/sergio/SPECTRA');
      fprintf(1,'%s \n',[homepath 'CO_COMMON/endpts.dat'])            
    end
    error('go back, cd to /home/sergio/SPECTRA and re run your code!');
  end

  fid = fopen(filename1,'w');
     
  a = strengthM*ones(1,3);
  b = [thebands' ll' hh'];

  fwrite(fid,a,'real*8');
  fwrite(fid,b,'real*8');

  fclose(fid);
end
